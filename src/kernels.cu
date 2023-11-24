#include <curand_kernel.h>
#include <cuda.h>
#include <math.h>

__device__ double abcd_prob(float a, float b, float c, float d, int logn, int num_ab){
    return pow(a+b, num_ab) * pow(c+d, logn-num_ab);
}

__device__ double ab_prob(float a, float b, int logn, int num_a){
    return pow(a, num_a) * pow(b, logn-num_a);
}

__device__ int bitCount(uint64_t vid){
    return __popcll(vid);
}

__device__ 