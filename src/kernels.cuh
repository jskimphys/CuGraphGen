#pragma once

#include <math.h>
#include <vector>
#include <thread>
#include <string>
#include <curand.h>
#include <cuda.h>

#include "SKGgenerator.cuh"
#include "fast_writer.h"

class CuWorker{
public:
    CuWorker(uint64_t seed);
    ~CuWorker();
    
    void process_workloads(std::vector<schedule_entry> workloads, std::string filename, size_t filesize, double a, double b, double c, double d);
    void update_random_arr(size_t num_32bit);
    size_t get_randomarr_bytesize(){
        return rarr_bytesize;
    }
    size_t get_edgearr_bytesize(){
        return earr_bytesize;
    }
private:
    size_t rarr_bytesize;
    size_t earr_bytesize;
    uint32_t* random_arr;
    uint64_t* edge_arr_device;
    curandGenerator_t gen;
    std::vector<cudaStream_t> streams;

    int hostMemIdx;
    std::vector<uint64_t*> edge_arr_host_list;
    fast_writer writer;
};
