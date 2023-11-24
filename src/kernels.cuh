#pragma once
#include "SKGgenerator.h"

uint64_t* cu_generate_edges(const schedule_entry& entry, uint16_t* randombits, double a, double b, double c, double d);
void process_workload(schedule_entry &entry, void* mem_start, uint64_t mem_size, uint64_t entry_seed, double a, double b, double c, double d);

class CUScheduler{
public:
    CUScheduler();
    ~CUScheduler();
    void CUScheduler::setup_cu_mem(uint64_t seed)
}