#include <curand_kernel.h>
#include <cuda.h>
#include <math.h>
#include <curand.h>
#include <vector>
#include <thread>

#include "kernels.cuh"

const int blockSize = 128;
const int STREAM_PER_WORKER = 16;

__global__ void generate_randombits_dst(uint64_t prefix, uint16_t bab, uint16_t dcd, int num_bits, uint16_t* random_array, uint64_t* output_array, uint64_t num_edges) {
    //inputs : after the bits in prefix, posterier bits are randomly generated
    //one_prob : probability of 1
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_edges) return;

    uint64_t rarr_idx = tid * num_bits;
    while(num_bits > 0) {
        if(output_array[2*tid] & (1 << num_bits)){//if the corresponding src_vid bit is 1
            if(random_array[rarr_idx++] < dcd){
                prefix = prefix | (1 << num_bits);
            }
        }
        else{
            if(random_array[rarr_idx++] < bab){
                prefix = prefix | (1 << num_bits);
            }
        }
        num_bits--;
    }
    output_array[2*tid+1] = prefix;
}

__global__ void generate_randombits_src(uint64_t prefix, uint16_t cdabcd, int num_bits, uint16_t* random_array, uint64_t* output_array, uint64_t num_edges) {
    //inputs : after the bits in prefix, posterier bits are randomly generated
    //cdabcd : probability of 1
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_edges) return;

    uint64_t rarr_idx = tid * num_bits;
    while(num_bits-- > 0) {
        if(random_array[rarr_idx++] < cdabcd){
            prefix = prefix | (1 << num_bits);
        }
    }
    output_array[2*tid] = prefix;
}

__global__ void fillWithStride2(uint64_t* data, uint64_t value, int size) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) {
        data[tid * 2] = value;
    }
}

CuWorker::CuWorker(uint64_t seed)
    : hostMemIdx(0), writer()
{
    size_t free_mem_size, total_mem_size;
    cudaMemGetInfo(&free_mem_size, &total_mem_size);
    
    //allocate memory for each worker
    size_t mem_per_worker = (size_t) (free_mem_size * 0.85);
    double rarr_earr_ratio = 4;

    const int host_mem_num = 2;//this is for writing parallelism
    earr_bytesize = mem_per_worker / (rarr_earr_ratio + 1);
    rarr_bytesize = mem_per_worker - earr_bytesize;

    for(int i = 0; i < STREAM_PER_WORKER; i++){
        cudaStream_t stream;
        cudaStreamCreate(&stream);
        streams.push_back(stream);
    }

    cudaMallocAsync((void**)&random_arr, rarr_bytesize + 8*1024, streams[0]);//extra 8KB for alignment(since wrongly aligned memory cause error)
    cudaMallocAsync((void**)&edge_arr_device,  earr_bytesize + 8*1024, streams[1]);
    
    
    edge_arr_host_list = std::vector<uint64_t*>(host_mem_num);
    for(int i=0; i< host_mem_num; i++){
        cudaMallocHost((void**)&(edge_arr_host_list[i]), earr_bytesize);
        if(cudaPeekAtLastError() != cudaSuccess){
            std::cerr << "Error in allocating memories" << std::endl;
            std::cerr << cudaGetErrorString(cudaPeekAtLastError()) << "at " << i << "th host memory allocation" << std::endl;
            exit(1);
        }
    }

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, seed);


    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "Error in allocating memories" << std::endl;
        std::cerr << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        exit(1);
    }
}

CuWorker::~CuWorker(){
    for(int i=0; i< edge_arr_host_list.size(); i++){
        cudaFreeHost(edge_arr_host_list[i]);
    }
    cudaFree(random_arr);
    cudaFree(edge_arr_device);
    curandDestroyGenerator(gen);
    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "Error in freeing memory or destroying rand generator" << std::endl;
        std::cerr << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        exit(1);
    }
    for(int i = 0; i < STREAM_PER_WORKER; i++){
        cudaStreamDestroy(streams[i]);
    }
}

void CuWorker::update_random_arr(size_t num_32bits){
    curandGenerate(gen, random_arr, num_32bits);
    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "Error in generating random bits" << std::endl;
        std::cerr << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        exit(1);
    }
}

void CuWorker::process_workloads(std::vector<schedule_entry> workloads, std::string filename, size_t filesize, double a, double b, double c, double d){

    size_t total_rbits = 0;
    size_t total_edges = 0;
    //caculate the total workload size
    for(auto entry : workloads){
        if(entry.t == schedule_entry::type::along_src_vid){
            total_rbits += entry.num_edge * (2*entry.log_n - entry.log_prefixlen);
        }
        else{
            total_rbits += entry.num_edge * (entry.log_n - entry.log_prefixlen);
        }
        total_edges += entry.num_edge;
    }
    update_random_arr(total_rbits / 2);//require 2 bytes of random value to generate one random bit in edge

    // maps memory for each workload
    uint64_t** edge_ptrs = new uint64_t*[workloads.size()];
    uint16_t** randombits_ptrs = new uint16_t*[workloads.size()];
    edge_ptrs[0] = (uint64_t*)edge_arr_device;
    randombits_ptrs[0] = (uint16_t*)random_arr;
    for(int i = 1; i < workloads.size(); i++){
        edge_ptrs[i] = edge_ptrs[i-1] + workloads[i-1].num_edge * 2;
        edge_ptrs[i] = edge_ptrs[i] + 16 - workloads[i-1].num_edge % 16;//align to 16 byest

        size_t randombits_needed = 0;
        if(workloads[i].t == schedule_entry::type::along_src_vid){
            //need to generate random bits for src_vid
            randombits_needed = workloads[i-1].num_edge * (workloads[i-1].log_n - workloads[i-1].log_prefixlen);
            randombits_ptrs[i] = randombits_ptrs[i-1] + randombits_needed;
            randombits_ptrs[i] = randombits_ptrs[i] + 16 - randombits_needed % 16;//align to 16 byest
            
            //need to generate random bits for dst_vid
            randombits_needed = workloads[i-1].num_edge * workloads[i-1].log_n;
            randombits_ptrs[i] = randombits_ptrs[i] + randombits_needed;
            randombits_ptrs[i] = randombits_ptrs[i] + 16 - randombits_needed % 16;//align to 16 byest
        }
        else{
            randombits_needed = workloads[i-1].num_edge * (workloads[i-1].log_n - workloads[i-1].log_prefixlen);
            randombits_ptrs[i] = randombits_ptrs[i-1] + randombits_needed;
            randombits_ptrs[i] = randombits_ptrs[i] + 16 - randombits_needed % 16;//align to 16 byest
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < workloads.size(); i++){

        schedule_entry& entry = workloads[i];
        if(entry.t == schedule_entry::type::along_dst_vid){
            //fill the src_vid in the edgelist
            int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
            uint16_t bab = (uint16_t) (b/(a+b) * (1 << 16));
            uint16_t dcd = (uint16_t) (d/(c+d) * (1 << 16));

            fillWithStride2<<<gridSize, blockSize, 0, streams[i % streams.size()]>>>((uint64_t*)edge_ptrs[i], entry.src_vid_start, entry.num_edge);
            generate_randombits_dst<<<gridSize, blockSize, 0, streams[i % streams.size()]>>>(entry.dst_vid_start, bab, dcd, entry.log_n - entry.log_prefixlen, randombits_ptrs[i], edge_ptrs[i], entry.num_edge);
        }
        else{
            int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
            uint16_t cdabcd = (uint16_t) ((c+d) * (1 << 16));
            generate_randombits_src<<<gridSize, blockSize, 0, streams[i % streams.size()]>>>(entry.src_vid_start, cdabcd, entry.log_n - entry.log_prefixlen, randombits_ptrs[i], edge_ptrs[i], entry.num_edge);

            
            uint16_t bab = (uint16_t) (b/(a+b) * (1 << 16));
            uint16_t dcd = (uint16_t) (d/(c+d) * (1 << 16));
            int randombits_used = entry.num_edge * (entry.log_n - entry.log_prefixlen);
            randombits_used = randombits_used + 16 - randombits_used % 16;//align to 16 byest
            generate_randombits_dst<<<gridSize, blockSize, 0, streams[i % streams.size()]>>>(entry.dst_vid_start, bab, dcd, entry.log_n, (uint16_t*)(randombits_ptrs[i] + randombits_used), (uint64_t*)edge_ptrs[i], entry.num_edge);
        }
    }


    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "Error processing workloads" << std::endl;
        std::cerr << cudaGetErrorString(cudaPeekAtLastError()) << std::endl;
        exit(1);
    }

    CUdeviceptr dvptr = CUdeviceptr(edge_arr_device);
    cuMemcpyDtoH(edge_arr_host_list[hostMemIdx], dvptr, total_edges * 16);
    writer.write_async(filename, (char*) edge_arr_host_list[hostMemIdx], total_edges * 16);
    hostMemIdx = (hostMemIdx + 1) % edge_arr_host_list.size();
}

// CuScheduler::CuScheduler(int num_workers){
//     //get cuda device properties
//     size_t free_mem_size, total_mem_size;
//     cudaMemGetInfo(&free_mem_size, &total_mem_size);
    
//     //allocate memory for each worker
//     size_t mem_per_worker = (size_t) (free_mem_size * 0.85 / num_workers);
//     double rarr_earr_ratio = 4;

//     earr_bytesize = mem_per_worker / (rarr_earr_ratio + 1);
//     rarr_bytesize = mem_per_worker - earr_bytesize;

//     this->num_workers = num_workers;
//     workers = new CuWorker*[num_workers];
//     for(int i = 0; i < num_workers; i++){
//         worker_status.push_back(false);
//         worker_threads.push_back(thread());
//         workers[i] = new CuWorker(i, rarr_bytesize, earr_bytesize, NUM_WRITER);
//     }
// }

// CuScheduler::~CuScheduler(){
//     //free memory for each worker
//     for(int i = 0; i < num_workers; i++){
//         worker_threads[i].join();
//     }
//     for(int i = 0; i < num_workers; i++){
//         delete workers[i];
//     }
// }

// void CuScheduler::deliver_workloads(std::vector<schedule_entry> workloads, std::string filename, size_t filesize, double a, double b, double c, double d){
//     //divide workloads to each worker
//     int idle_worker = -1;
//     for(int i=0; i < num_workers; i++){
//         if(worker_status[i] == false){
//             idle_worker = i;
//             break;
//         }
//     }
//     if(idle_worker == -1){
//         while(true){//spin wait for idle worker
//             for(int i=0; i<num_workers; i++){
//                 if(worker_threads[i].joinable()){
//                     worker_threads[i].join();
//                     worker_status[i] = false;
//                     idle_worker = i;
//                     break;
//                 }
//             }
//             if(idle_worker != -1) break;
//         }
//     }

//     worker_threads[idle_worker] = thread(&CuWorker::process_workloads, workers[idle_worker], workloads, filename, filesize, a, b, c, d);
//     worker_status[idle_worker] = true;
// }

// size_t CuScheduler::get_randomarr_bytesize(){
//     return rarr_bytesize;
// }
// size_t CuScheduler::get_edgearr_bytesize(){
//     return earr_bytesize;
// }