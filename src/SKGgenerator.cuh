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
    SKGgenerator(double a, double b, double c, double d, int log_n, int edge_ratio, uint64_t seed, std::string dir, uint64_t workload_byte_limit) : a(a), b(b), c(c), d(d), log_n(log_n), seed(seed), dir(dir), workload_byte_limit(workload_byte_limit){
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
        n_edge = n_vertex * edge_ratio;
        if (n_edge > n_vertex * n_vertex){
            std::cout << "This algorithm only works well for spase graph. since Each edge is generated in independently, |E|~|V|^2 case will give many overlapping edges" << std::endl;
            exit(0);
        }
    }
    void generate();
    void divide_workloads(); //generates workload by devide original workload
    void divide_workloads_naive();
    void print_workload_summary(){
        std::cout << "-----workload generation summary-----" << std::endl;
        std::cout << "number of workloads: " << workloads.size() << std::endl;
        uint64_t total_edge = 0;
        for (int i = 0; i < workloads.size(); i++){
            total_edge += workloads[i].num_edge;
        }

        std::cout << "total vertex     : " << n_vertex << std::endl;
        std::cout << "target total edge: " << n_edge << std::endl;
        std::cout << "real total edge  : " << total_edge << std::endl;
        std::cout << "-------------------------------------" << std::endl;
        n_edge = total_edge;//this might be different from the original n_edge since we flunctuate the number of edges in each workload
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
    uint64_t n_vertex;
    uint64_t n_edge;

    std::string dir;

    uint64_t seed;
    uint64_t workload_byte_limit;
    std::vector<schedule_entry> workloads;

    uint64_t filesize_limit = 1ULL << 30; // 1GB
};

