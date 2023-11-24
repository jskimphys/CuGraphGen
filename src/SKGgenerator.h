#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <bitset>

struct schedule_entry {
    enum class type{
        along_src_vid,
        along_dst_vid,
    };
    type t;
    uint64_t src_vid_start;
    uint64_t src_vid_end;
    uint64_t dst_vid_start;
    uint64_t dst_vid_end;
    uint64_t num_edge;
    int log_n;
    int log_prefixlen;
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

        n_vertex = (1ULL << log_n);
        n_edge = (1ULL << log_edge);
        if (n_edge > n_vertex * n_vertex){
            std::cout << "Too many edges, thus limit it to n_vertex square" << std::endl;
            n_edge = n_vertex * n_vertex;
        }
    }
    void generate();
    //generates workload by devide original workload
    void schedule();
    void print_workload_summary(){
        std::cout << "-----workload generation summary-----" << std::endl;
        std::cout << "number of workloads: " << workloads.size() << std::endl;
        uint64_t total_edge = 0;
        for (int i = 0; i < workloads.size(); i++){
            total_edge += workloads[i].num_edge;
        }

        std::cout << "total vertex: " << n_vertex << std::endl;
        std::cout << "target total edge: " << n_edge << std::endl;
        std::cout << "real total edge  : " << total_edge << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        n_edge = total_edge;//this might be different from the original n_edge since we flunctuate the number of edges in each workload
    }
    void print_workload(){//this is for debug
        for (int i = 0; i < workloads.size(); i++){
            std::cout << "workload " << i << std::endl;
            //print in binary
            std::cout << "src_vid_start: " << std::bitset<64>(workloads[i].src_vid_start) << std::endl;
            std::cout << "src_vid_end  : " << std::bitset<64>(workloads[i].src_vid_end) << std::endl;
            std::cout << "dst_vid_start: " << std::bitset<64>(workloads[i].dst_vid_start) << std::endl;
            std::cout << "dst_vid_end  : " << std::bitset<64>(workloads[i].dst_vid_end) << std::endl;
            std::cout << "num_edge     : " << std::bitset<64>(workloads[i].num_edge) << std::endl;
        }
    }
private:
    uint64_t workload_size_calc_src_dim(int log_n, uint64_t src_vid_start, uint64_t src_vid_end);
    uint64_t workload_size_calc_dst_dim(int log_n, uint64_t dst_vid_start, uint64_t dst_vid_end, uint64_t num_edge_in_src_dim);
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

    uint64_t filesize_limit = 1ULL << 30; // 1GB
};

int get_lsb_loc(uint64_t n);
int get_msb_loc(uint64_t n);
int count_bits(uint64_t n);
uint64_t flunctuate(uint64_t n, double prob);