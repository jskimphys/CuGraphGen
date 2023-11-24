#include <cmath>
#include <random>
#include <cassert>
#include <bitset>

#include "SKGgenerator.h"

using namespace std;

//schedule workloads by divide original workload
//devision is done along the src vertex(conceptually row of the adjecency matrix)
//so each workload is a submatrix with a comparably small number of rows and all columns
//each division is only happen at vid = 2^k to make it easier to consume the workload
void SKGgenerator::schedule(){
    if(a+b < c+d){
        std::swap(a, c);
        std::swap(b, d);
    }

    uint64_t src_vid_start = 0;
    uint64_t src_vid_end = n_vertex;
    uint64_t max_edge_per_workload = workload_size_limit / 8;// 8 is the size of a edge in bytes

    while(src_vid_start < n_vertex){
        int log_postfix = get_lsb_loc(src_vid_end - src_vid_start);
        int bitcount = count_bits(src_vid_start >> log_postfix);
        double prob = pow(c+d, bitcount) * pow(a+b, log_n - log_postfix - bitcount);
        uint64_t num_edge_in_workload =  static_cast<uint64_t>(n_edge * prob);

        //print in binary
        cout << "src_vid_start: " << bitset<64>(src_vid_start) << endl;
        cout << "src_vid_end  : " << bitset<64>(src_vid_end) << endl;
        //not now
        cout << "prob: " << prob << endl;
        cout << "num_edge_in_workload: " << num_edge_in_workload << endl;

        if (num_edge_in_workload > max_edge_per_workload){
            assert(src_vid_end - src_vid_start != 1);
            src_vid_end = src_vid_end - (src_vid_end - src_vid_start) / 2;
            continue;
        }

        schedule_entry entry{
            src_vid_start,
            src_vid_end,
            0,
            n_vertex,
            num_edge_in_workload
        };
        
        workloads.push_back(entry);
        src_vid_start = src_vid_end;
        src_vid_end = src_vid_start +  1LL << get_lsb_loc(src_vid_start);
    }
    return;
}

int get_lsb_loc(uint64_t n){
    int count = -1;
    while(n){
        count++;
        n >>= 1;
    }
    return count;
}
int get_msb_loc(uint64_t n){
    int count = 0;
    while(n){
        count++;
        n >>= 1;
    }
    return count;
}

int count_bits(uint64_t n){
    int count = 0;
    while(n){
        if(n & 1){
            count++;
        }
        n >>= 1;
    }
    return count;
}


uint64_t flunctuate(uint64_t n_vertex, double prob){
    double stdev = sqrt(n_vertex * prob * (1 - prob));
    double mean = n_vertex * prob;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(mean, stdev);

    uint64_t result = static_cast<uint64_t>(d(gen));
    return result;
}