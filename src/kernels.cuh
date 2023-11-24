#pragma once
#include <math.h>
#include <vector>
#include "SKGgenerator.h"

void deliver_workloads(std::vector<schedule_entry> workloads, void* mem_start, double a, double b, double c, double d);
void start_scheduler(uint64_t seed);
void terminate_scheduler();

size_t get_randomarr_size();
size_t get_edgearr_size();