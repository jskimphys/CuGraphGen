#include <curand_kernel.h>
#include <cuda.h>
#include <math.h>
#include <curand.h>
#include <vector>

#include "kernels.cuh"


class CuScheduler{
public:
    CuScheduler(){
        random_arr = nullptr;
        edge_arr = nullptr;
        random_arr_size = 0;
        edge_arr_size = 0;
    }
    void setup_cu_mem(uint64_t seed);
    void free_cu_mem();
    void process_workloads(std::vector<schedule_entry> workloads, void* mem_start, double a, double b, double c, double d);
    void update_random_arr(size_t update_size);
    size_t get_randomarr_size(){
        return random_arr_size;
    }
    size_t get_edgearr_size(){
        return edge_arr_size;
    }
private:
    size_t random_arr_size;
    size_t edge_arr_size;
    uint32_t* random_arr;
    uint64_t* edge_arr;
    curandGenerator_t gen;
} cu_scheduler;

__global__ void generate_randombits_dst(uint64_t prefix, uint16_t one_prob, uint64_t num_bits, uint16_t* random_array, uint64_t* output_array, uint64_t num_edges) {
    //inputs : after the bits in prefix, posterier bits are randomly generated
    //one_prob : probability of 1
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_edges) return;

    uint64_t rarr_idx = tid * num_bits;
    while(num_bits-- > 0) {
        if(random_array[rarr_idx++] < one_prob){
            prefix = prefix | (1 << num_bits);
        }
    }
    output_array[2*tid+1] = prefix;
}

__global__ void generate_randombits_src(uint64_t prefix, uint16_t one_prob, uint64_t num_bits, uint16_t* random_array, uint64_t* output_array, uint64_t num_edges) {
    //inputs : after the bits in prefix, posterier bits are randomly generated
    //one_prob : probability of 1
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_edges) return;

    uint64_t rarr_idx = tid * num_bits;
    while(num_bits-- > 0) {
        if(random_array[rarr_idx++] < one_prob){
            prefix = prefix | (1 << num_bits);
        }
    }
    output_array[2*tid] = prefix;
}

__global__ void fillWithStride2(uint64_t* data, uint64_t value, int size) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) {
        // Calculate the offset with a stride of 2
        int offset = tid * 2;
        data[offset] = value;
    }
}

void CuScheduler::setup_cu_mem(uint64_t seed){
    size_t freeBytes, totalBytes;
    cudaMemGetInfo(&freeBytes, &totalBytes);

    std::cout << "Free cuda memory: " << freeBytes / (1024 * 1024) << " MB" << std::endl;
    std::cout << "Total cuda memory: " << totalBytes / (1024 * 1024) << " MB" << std::endl;

    size_t memory_to_use = (size_t) (freeBytes * 0.85);
    //divide memory into random bits and edges
    //edge_size = 16*n_edges
    //random_bits_size = 2*random_bits_per_edge * n_edges ~ 2*40*n_edges
    int rbits_edges_ratio = 5;
    edge_arr_size = memory_to_use / (rbits_edges_ratio + 1);
    random_arr_size = memory_to_use - edge_arr_size;

    cudaMalloc((void**)&random_arr, random_arr_size);
    cudaMalloc((void**)&edge_arr, edge_arr_size);
    if(random_arr == nullptr || edge_arr == nullptr){
        std::cerr << "Error in allocating memory for random bits or edges" << std::endl;
        exit(1);
    }

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, seed);
}

void CuScheduler::free_cu_mem(){
    cudaFree(random_arr);
    cudaFree(edge_arr);
    curandDestroyGenerator(gen);
}

void CuScheduler::update_random_arr(size_t update_size){
    curandGenerate(gen, random_arr, random_arr_size / sizeof(uint32_t));
    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "Error in generating random bits" << std::endl;
        exit(1);
    }
}

void CuScheduler::process_workloads(std::vector<schedule_entry> workloads, void* mem_start, double a, double b, double c, double d){
    size_t total_bits_needed = 0;
    size_t total_edges = 0;
    for(auto entry : workloads){
        if(entry.t == schedule_entry::type::along_src_vid){
            total_bits_needed += entry.num_edge * (entry.log_n - entry.log_prefixlen);
        }
        else{
            total_bits_needed += entry.num_edge * (2*entry.log_n - entry.log_prefixlen);
        }
        total_edges += entry.num_edge;
    }
    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "somthing got wrong before generating random bits" << std::endl;
        exit(1);
    }
    update_random_arr(total_bits_needed >> 1 + 1);

    char** edge_ptrs = new char*[workloads.size()];
    char** randombits_ptrs = new char*[workloads.size()];
    edge_ptrs[0] = (char*)edge_arr;
    randombits_ptrs[0] = (char*)random_arr;
    for(int i = 1; i < workloads.size(); i++){
        edge_ptrs[i] = edge_ptrs[i-1] + workloads[i-1].num_edge * 16;
        size_t randombits_needed = (workloads[i-1].t == schedule_entry::type::along_src_vid) ? workloads[i-1].num_edge * (workloads[i-1].log_n - workloads[i-1].log_prefixlen) : workloads[i-1].num_edge * (2*workloads[i-1].log_n - workloads[i-1].log_prefixlen);
        randombits_ptrs[i] = randombits_ptrs[i-1] + randombits_needed;
        randombits_ptrs[i] = randombits_ptrs[i] + 16 - randombits_needed % 16;
    }

    #pragma omp parallel for num_threads(workloads.size() / 64)
    for(int i = 0; i < workloads.size(); i++){
        schedule_entry& entry = workloads[i];
        if(entry.t == schedule_entry::type::along_dst_vid){
            //fill the src_vid in the edgelist
            int blockSize = 256;
            int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
            uint16_t one_prob = (uint16_t) ((b+d) * (1 << 16));//probability of 1 quantized to 16 bits
            fillWithStride2<<<gridSize, blockSize>>>((uint64_t*)edge_ptrs[i], entry.src_vid_start, entry.num_edge);

            generate_randombits_dst<<<gridSize, blockSize>>>(entry.dst_vid_start, one_prob, entry.log_n - entry.log_prefixlen, (uint16_t*)randombits_ptrs[i], (uint64_t*)edge_ptrs[i], entry.num_edge);
        }
        else{
            int blockSize = 256;
            int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
            uint16_t one_prob = (uint16_t) ((c+d) * (1 << 16));
            generate_randombits_src<<<gridSize, blockSize>>>(entry.src_vid_start, one_prob, entry.log_n - entry.log_prefixlen, (uint16_t*)randombits_ptrs[i], (uint64_t*)edge_ptrs[i], entry.num_edge);

            one_prob = (uint16_t) ((b+d) * (1 << 16));
            int randombits_needed = entry.num_edge * entry.log_n  + 16 - entry.num_edge * entry.log_n % 16;
            generate_randombits_dst<<<gridSize, blockSize>>>(entry.dst_vid_start, one_prob, entry.log_n, (uint16_t*)randombits_ptrs[i] + randombits_needed, (uint64_t*)edge_ptrs[i], entry.num_edge);
        }
    }

    if(cudaPeekAtLastError() != cudaSuccess){
        std::cerr << "somthing got wrong before memcpy" << std::endl;
        exit(1);
    }

    cudaDeviceSynchronize();
    CUdeviceptr cu_edge_ptr = CUdeviceptr(edge_arr);
    cuMemcpyDtoH(mem_start, cu_edge_ptr, total_edges * 16);
}


void deliver_workloads(std::vector<schedule_entry> workloads, void* mem_start, double a, double b, double c, double d){
    cu_scheduler.process_workloads(workloads, mem_start, a, b, c, d);
}

void start_scheduler(uint64_t seed){
    cu_scheduler.setup_cu_mem(seed);
}

void terminate_scheduler(){
    cu_scheduler.free_cu_mem();
}

size_t get_randomarr_size(){
    return cu_scheduler.get_randomarr_size();
}
size_t get_edgearr_size(){
    return cu_scheduler.get_edgearr_size();
}