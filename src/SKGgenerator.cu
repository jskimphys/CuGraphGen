#include <cmath>
#include <random>
#include <cassert>
#include <bitset>
#include <filesystem>

#include "kernels.cuh"
#include "SKGgenerator.cuh"

using namespace std;

const int EDGE_BYTE = 16;

int get_lsb_loc(uint64_t n);
int get_msb_loc(uint64_t n);
int count_bits(uint64_t n);
int64_t flunctuate(uint64_t n, double prob, mt19937& gen);


//schedule workloads by divide original workload
//devision is done along the src vertex(conceptually row of the adjecency matrix)
//so each workload is a submatrix with a comparably small number of rows and all columns
//each division is only happen at vid = 2^k to make it easier to consume the workload
void SKGgenerator::divide_workloads(){
    if(a+b < c+d){//this part is just to make workload division more straightforward
        std::swap(a, c);
        std::swap(b, d);
    }
    // assert(a+c >= b+d);

    uint64_t src_vid_start = 0;
    uint64_t src_vid_end = n_vertex;
    uint64_t max_edge_per_workload = (uint64_t)(workload_byte_limit / EDGE_BYTE * 0.85);
    double edge_prob = 1;
    
    while(src_vid_start < n_vertex){
        if (n_edge * edge_prob > max_edge_per_workload){
            if(src_vid_end - src_vid_start == 1){
                //if src vid is all divided, but desirable size is not reached
                //now divide in dst_vid dimension
                uint64_t dst_vid_start = 0;
                uint64_t dst_vid_end = n_vertex;
                while(dst_vid_start < n_vertex){
                    if (n_edge * edge_prob > max_edge_per_workload){
                        dst_vid_end = dst_vid_end - (dst_vid_end - dst_vid_start) / 2;
                        edge_prob *= (a+c);
                        continue;
                    }

                    schedule_entry entry{
                        .t = schedule_entry::type::along_dst_vid,
                        .src_vid_start = src_vid_start,
                        .src_vid_end = src_vid_end,
                        .dst_vid_start = dst_vid_start,
                        .dst_vid_end = dst_vid_end,
                        .num_edge = (uint64_t)round(n_edge * edge_prob),
                        .log_n = log_n,
                        .log_prefixlen = get_lsb_loc(dst_vid_end - dst_vid_start)
                    };
                    
                    workloads.push_back(entry);
                    dst_vid_start = dst_vid_end;
                    dst_vid_end = dst_vid_start +  1ULL << get_lsb_loc(dst_vid_start);
                }
                src_vid_start = src_vid_end;
                src_vid_end = src_vid_start +  1ULL << get_lsb_loc(src_vid_start);
                continue;
            }
            src_vid_end = src_vid_end - (src_vid_end - src_vid_start) / 2;
            edge_prob *= (a+b);
            continue;
        }

        schedule_entry entry{
            .t = schedule_entry::type::along_src_vid,
            .src_vid_start = src_vid_start,
            .src_vid_end = src_vid_end,
            .dst_vid_start = 0,
            .dst_vid_end = n_vertex,
            .num_edge = (uint64_t)round(n_edge * edge_prob),
            .log_n = log_n,
            .log_prefixlen = get_lsb_loc(src_vid_end - src_vid_start)
        };
        
        workloads.push_back(entry);
        if(get_lsb_loc(src_vid_start) > get_lsb_loc(src_vid_end)){
            //ex) current
            //src_vid_start = 001000
            //src_vid_end   = 001010
            src_vid_start = src_vid_end;
            src_vid_end = src_vid_end + (src_vid_end - src_vid_start);
            edge_prob *= (c+d)/(a+b);
        }else{
            //ex) current
            //src_vid_start = 00111110
            //src_vid_end   = 01000000
            src_vid_start = src_vid_end;
            src_vid_end = src_vid_end + (src_vid_end - src_vid_start);
            edge_prob *= (a+b)/(c+d);
            //current state
            //src_vid_start = 01000000
            //src_vid_end   = 01000010
            while(get_lsb_loc(src_vid_start) > get_lsb_loc(src_vid_end)){
                src_vid_start = src_vid_end;
                src_vid_end = src_vid_end + (src_vid_end - src_vid_start);
                edge_prob /= (c+d);
            }
        }
    }

    //since probability is not exactly the same as the number of edges
    //we flunctuate the number of edges in each workload
    int wksize = workloads.size();
    random_device rd;
    mt19937 gen(rd());
    gen.seed(seed);

    for(int i = 0; i < wksize; i++){
        schedule_entry& entry = workloads[i];
        entry.num_edge = max(0L, flunctuate(n_edge, (double) entry.num_edge / (double) n_edge, gen));
    }
    return;
}

void SKGgenerator::generate(){
    if(filesystem::is_directory(dir)){
        std::cout << "Directory already exists, deleting it" << std::endl;
        filesystem::remove_all(dir);
    }
    filesystem::create_directory(dir);
    
    if(filesystem::space(dir).available < n_edge * 18){
        std::cout << "available storage : " << filesystem::space(dir).available << std::endl;
        std::cout << "required storage  : " << n_edge * 18 << std::endl;
        std::cout << "Not enough storage, exiting" << std::endl;
        exit(0);
    }


    CuWorker cuworker(seed);
    size_t earr_bytesize = cuworker.get_edgearr_bytesize();
    size_t rarr_bytesize = cuworker.get_randomarr_bytesize();


    struct file_info{
        int file_id;
        uint64_t size;
        vector<schedule_entry> mappings;
    };
    vector<file_info> file_infos;
    

    //map workload entries to files
    size_t infile_address = 0;
    size_t needed_random_arr_size = 0;
    int file_id = 0;
    file_infos.push_back(file_info{file_id, 0});
    for(auto& entry : workloads){
        size_t entry_random_arr_size = 0;
        if(entry.t == schedule_entry::type::along_src_vid){
            entry_random_arr_size = entry.num_edge * (2*entry.log_n - entry.log_prefixlen) * 2;
        }else{
            entry_random_arr_size = entry.num_edge * (entry.log_n - entry.log_prefixlen) * 2;
        }

        if(infile_address + entry.num_edge * 16 > earr_bytesize || needed_random_arr_size + entry_random_arr_size > rarr_bytesize){
            file_infos[file_id].size = infile_address;
            infile_address = entry.num_edge * 16;
            needed_random_arr_size = entry_random_arr_size;

            file_id++;
            file_infos.push_back(file_info{file_id, 0});
        }
        file_infos[file_id].mappings.push_back(entry);
        infile_address += entry.num_edge * 16;
        needed_random_arr_size += entry_random_arr_size;
    }
    
    file_infos[file_id].size = infile_address;

    //get memory and start to process each workload

    for(auto& file_info : file_infos){
        cuworker.process_workloads(file_info.mappings, dir + "/edgelist8B_" + to_string(file_info.file_id) + ".part", file_info.size, a, b, c, d);
    }
}

uint64_t SKGgenerator::workload_size_calc_src_dim(int log_n, uint64_t src_vid_start, uint64_t src_vid_end){
    int log_postfix = get_lsb_loc(src_vid_end - src_vid_start);
    int bitcount = count_bits(src_vid_start >> log_postfix);
    double prob = pow(c+d, bitcount) * pow(a+b, log_n - log_postfix - bitcount);
    //assert(prob <= 1);
    uint64_t num_edge_in_workload =  static_cast<uint64_t>(round(n_edge * prob));
    return num_edge_in_workload;
}

uint64_t SKGgenerator::workload_size_calc_dst_dim(int log_n, uint64_t dst_vid_start, uint64_t dst_vid_end, uint64_t num_edge_in_src_dim){
    int log_postfix = get_lsb_loc(dst_vid_end - dst_vid_start);
    int bitcount = count_bits(dst_vid_start >> log_postfix);
    double prob = pow(b+d, bitcount) * pow(a+c, log_n - log_postfix - bitcount);
    uint64_t num_edge_in_workload =  static_cast<uint64_t>(round(num_edge_in_src_dim * prob));
    return num_edge_in_workload;
}

int get_lsb_loc(uint64_t n){
    if(n == 0){
        return 64;
    }
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
    if (n == 0){
        return 0;
    }
    int count = 0;
    while(n){
        if(n & 1){
            count++;
        }
        n >>= 1;
    }
    return count;
}


inline int64_t flunctuate(uint64_t n, double prob, mt19937& gen){
    double stdev = sqrt(n * prob * (1 - prob));
    double mean = n * prob;

    normal_distribution<> d(mean, stdev);

    int64_t result = static_cast<int64_t>(d(gen));
    return result;
}