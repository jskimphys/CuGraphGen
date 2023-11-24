#include <curand_kernel.h>
#include <cuda.h>
#include <math.h>
#include <curand.h>
#include <vector>

#include "kernels.cuh"

using namespace std;

uint* random_arr;
uint64_t* edges_in_cuda;
curandGenerator_t gen;

__global__ void generate_randombits_dst(uint64_t prefix, uint16_t one_prob, uint64_t num_bits, uint16_t* random_array, uint64_t* output_array, uint64_t num_edges) {
    //inputs : after the bits in prefix, posterier bits are randomly generated
    //one_prob : probability of 1
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if(tid >= num_edges) return;

    uint64_t rarr_idx = tid * num_bits;
    uint64_t post_bits = 0;
    while(num_bits-- > 0) {
        if(random_array[rarr_idx] < one_prob){
            prefix = prefix | (1 << post_bits);
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
    uint64_t post_bits = 0;
    while(num_bits-- > 0) {
        if(random_array[rarr_idx] < one_prob){
            prefix = prefix | (1 << post_bits);
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


uint64_t* cu_generate_edges(const schedule_entry& entry, uint16_t* randombits, double a, double b, double c, double d){
    uint64_t* edges = nullptr;
    auto err = cudaMalloc((void**)&edges, entry.num_edge * 2*sizeof(uint64_t));
    if(err != cudaSuccess) {
        std::cerr << "Error in allocating memory for edges" << std::endl;
        exit(1);
    }

    if(entry.t == schedule_entry::type::along_dst_vid){
        //fill the src_vid in the edgelist
        int blockSize = 1024;
        int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
        fillWithStride2<<<gridSize, blockSize>>>(edges, entry.src_vid_start, entry.num_edge);

        //fill the dst_vid in the edgelist
        blockSize = 256;
        gridSize = (entry.num_edge + blockSize - 1) / blockSize;
        uint16_t one_prob = (uint16_t) ((b+d) * (1 << 16));//probability of 1 quantized to 16 bits
        generate_randombits_dst<<<gridSize, blockSize>>>(entry.dst_vid_start, one_prob, entry.log_n - entry.log_prefixlen, randombits, edges, entry.num_edge);
    }
    else{
        //fill the src_vid in the edgelist
        int blockSize = 256;
        int gridSize = (entry.num_edge + blockSize - 1) / blockSize;
        uint16_t one_prob = (uint16_t) ((c+d) * (1 << 16));//probability of 1 quantized to 16 bits
        generate_randombits_src<<<gridSize, blockSize>>>(entry.src_vid_start, one_prob, entry.log_n - entry.log_prefixlen, randombits, edges, entry.num_edge);

        //fill the dst_vid in the edgelist
        one_prob = (uint16_t) ((b+d) * (1 << 16));//probability of 1 quantized to 16 bits
        generate_randombits_dst<<<gridSize, blockSize>>>(entry.dst_vid_start, one_prob, entry.log_n - entry.log_prefixlen, randombits, edges, entry.num_edge);
    }   

    return edges;
}

void CUScheduler::setup_cu_mem(uint64_t seed){
    size_t freeBytes, totalBytes;
    cudaMemGetInfo(&freeBytes, &totalBytes);

    cout << "Free cuda memory: " << freeBytes / (1024 * 1024) << " MB" << endl;
    cout << "Total cuda memory: " << totalBytes / (1024 * 1024) << " MB" << endl;

    size_t memory_to_use = (size_t) (freeBytes * 0.85);
    //divide memory into random bits and edges
    //edge_size = 16*n_edges
    //random_bits_size = 2*random_bits_per_edge * n_edges ~ 2*24*n_edges
    int rbits_edges_ratio = 3;
    size_t random_bits_size = memory_to_use / (rbits_edges_ratio + 1);
    size_t edge_size = memory_to_use - random_bits_size;

    cudaMalloc((void**)&random_arr, random_bits_size);
    cudaMalloc((void**)&edges_in_cuda, edge_size);

    curandCreateGenerator(&gen, CURAND_RNG_PSEUDO_DEFAULT);
    curandSetPseudoRandomGeneratorSeed(gen, seed);
}

void free_cu_mem(){
    cudaFree(random_arr);
    cudaFree(edges_in_cuda);
    curandDestroyGenerator(gen);
}

void process_workload(schedule_entry &entry, void* mem_start, uint64_t mem_size, uint64_t entry_seed, double a, double b, double c, double d){
    uint32_t random_uint32_needed;
    if(entry.t == schedule_entry::type::along_src_vid){
        //use one 32bit random number to generate 2 16bit random numbers
        random_uint32_needed = entry.num_edge * (2*entry.log_n - entry.log_prefixlen) / 2;
    }
    else{
        //use one 32bit random number to generate 2 16bit random numbers
        random_uint32_needed = entry.num_edge * (entry.log_n - entry.log_prefixlen) / 2;
    }
        
    curandGenerator_t gen;
    curandGenerate(gen, random_arr, random_uint32_needed);

    edges_in_cuda = cu_generate_edges(entry, (uint16_t*)random_arr, a, b, c, d);
    CUdeviceptr edges_in_cuda_ptr = CUdeviceptr(edges_in_cuda);
    cuMemcpyDtoH(mem_start, edges_in_cuda_ptr, entry.num_edge * 16);
}