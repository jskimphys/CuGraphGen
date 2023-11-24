#pragma once
#include "graph.h"
#include <string>
#include <vector>
#include <iostream>


struct schedule_entry {
    uint64_t src_vid_start;
    uint64_t src_vid_end;
    uint64_t dst_vid_start;
    uint64_t dst_vid_end;
    uint64_t num_edge;
};

class SKGgenerator {
public:
    SKGgenerator(double a, double b, double c, double d, int log_n, int log_edge, uint64_t seed, std::string dir, uint64_t workload_size_limit) : a(a), b(b), c(c), d(d), log_n(log_n), log_edge(log_edge), seed(seed), dir(dir), workload_size_limit(workload_size_limit){
        if (a +b < c + d){
            std::swap(a, c);
            std::swap(b, d);
        }
        double norm = a + b + c + d;
        a /= norm;
        b /= norm;
        c /= norm;
        d /= norm;

        n_vertex = (1LL << log_n);
        n_edge = (1LL << log_edge);
    }
    void generate();
    //generates workload by devide original workload
    void schedule();
    void print_schedule(){
        for (int i = 0; i < workloads.size(); i++){
            std::cout << "workload " << i << std::endl;
            std::cout << "src_vid_start: " << workloads[i].src_vid_start << std::endl;
            std::cout << "src_vid_end: " << workloads[i].src_vid_end << std::endl;
            std::cout << "dst_vid_start: " << workloads[i].dst_vid_start << std::endl;
            std::cout << "dst_vid_end: " << workloads[i].dst_vid_end << std::endl;
            std::cout << "num_edge: " << workloads[i].num_edge << std::endl;
        }
    }
    void process_workload(schedule_entry& entry);
private:
    double a;
    double b;
    double c;
    double d;
    int log_n;
    int log_edge;
    uint64_t n_vertex;
    uint64_t n_edge;

    std::string dir;

    uint64_t seed;
    uint64_t workload_size_limit;
    std::vector<schedule_entry> workloads;
};

int get_lsb_loc(uint64_t n);
int get_msb_loc(uint64_t n);
int count_bits(uint64_t n);
uint64_t flunctuate(uint64_t n, double prob);