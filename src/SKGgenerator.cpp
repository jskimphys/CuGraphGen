#include <cmath>
#include <random>
#include <cassert>
#include <bitset>
#include <filesystem>

#include "kernels.cuh"
#include "SKGgenerator.h"
#include "large_file.h"

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
        uint64_t num_edge_in_workload = workload_size_calc_src_dim(log_n, src_vid_start, src_vid_end);

        if (num_edge_in_workload > max_edge_per_workload){
            //cout << "workload size too large, divide it : " << num_edge_in_workload << endl;

            if(src_vid_end - src_vid_start == 1){
                //now divide in dst_vid dimension
                uint64_t dst_vid_start = 0;
                uint64_t dst_vid_end = n_vertex >> 1;
                while(dst_vid_start < n_vertex){
                    uint64_t num_edge_in_dst_workload = workload_size_calc_dst_dim(log_n, dst_vid_start, dst_vid_end, num_edge_in_workload);
                    if (num_edge_in_dst_workload > max_edge_per_workload){
                        dst_vid_end = dst_vid_end - (dst_vid_end - dst_vid_start) / 2;
                        continue;
                    }

                    schedule_entry entry{
                        .t = schedule_entry::type::along_dst_vid,
                        .src_vid_start = src_vid_start,
                        .src_vid_end = src_vid_end,
                        .dst_vid_start = dst_vid_start,
                        .dst_vid_end = dst_vid_end,
                        .num_edge = num_edge_in_dst_workload,
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
            continue;
        }

        schedule_entry entry{
            .t = schedule_entry::type::along_src_vid,
            .src_vid_start = src_vid_start,
            .src_vid_end = src_vid_end,
            .dst_vid_start = 0,
            .dst_vid_end = n_vertex,
            .num_edge = num_edge_in_workload,
            .log_n = log_n,
            .log_prefixlen = get_lsb_loc(src_vid_end - src_vid_start)
        };
        
        workloads.push_back(entry);
        src_vid_start = src_vid_end;
        src_vid_end = src_vid_start +  (1ULL << get_lsb_loc(src_vid_start));
    }

    //since probability is not exactly the same as the number of edges
    //we flunctuate the number of edges in each workload
    int wksize = workloads.size();
    #pragma omp parallel for 
    for(int i = 0; i < wksize; i++){
        schedule_entry& entry = workloads[i];
        entry.num_edge = entry.num_edge + flunctuate(n_vertex, (double) entry.num_edge / (double) n_edge);
    }
    return;
}

void SKGgenerator::generate(){
    if(filesystem::is_directory(dir)){
        std::cout << "Directory already exists, deleting it" << std::endl;
        filesystem::remove_all(dir);
    }
    filesystem::create_directory(dir);
    
    if(get_avail_storage(dir) < n_edge * 14){
        std::cout << "Not enough storage, exiting" << std::endl;
        exit(0);
    }

    start_scheduler(seed);
    size_t random_arr_size = get_randomarr_size();
    size_t edge_arr_size = get_edgearr_size();
    cout << "random arr size: " << random_arr_size << endl;
    cout << "edge arr size: " << edge_arr_size << endl;

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
        //cout << entry.num_edge*16 << endl;
        //cout << filesize_limit << endl;
        if(entry.t == schedule_entry::type::along_src_vid){
            needed_random_arr_size += entry.num_edge * (2*entry.log_n - entry.log_prefixlen) * 2;
        }else{
            needed_random_arr_size += entry.num_edge * (entry.log_n - entry.log_prefixlen) * 2;
        }
        // cout << "needed random arr size: " << needed_random_arr_size << endl;

        if(infile_address + entry.num_edge * 16 > edge_arr_size || needed_random_arr_size > random_arr_size){
            // if(needed_random_arr_size > random_arr_size){
            //     cout << "random arr size not enough" << endl;
            // }
            // else{
            //     cout << "edge arr size limit reached" << endl;
            // }
            file_infos[file_id].size = infile_address;
            infile_address = 0;
            file_id++;
            file_infos.push_back(file_info{file_id, 0});
            needed_random_arr_size = 0;
        }
        file_infos[file_id].mappings.push_back(entry);
        infile_address += entry.num_edge * 16;
    }
    
    file_infos[file_id].size = infile_address;

    //get memory and start to process each workload

    //print summary of mappings
    for(auto& file_info : file_infos){
        // for(auto& mapping : file_info.mappings){
        //     cout << "    workload " << mapping.entry.src_vid_start << " " << mapping.entry.src_vid_end << " " << mapping.entry.dst_vid_start << " " << mapping.entry.dst_vid_end << " " << mapping.entry.num_edge << endl;
        // }
    }

    for(auto& file_info : file_infos){
        char* mem = (char*)get_mmap_memory(dir + "/edgelist8B_" + to_string(file_info.file_id) + ".part", file_info.size);
            //void deliver_workloads(vector<schedule_entry> workloads, void* mem_start, double a, double b, double c, double d);
        deliver_workloads(file_info.mappings, mem, a, b, c, d);
        free_mmap_memory(mem, file_info.size);
    }
    terminate_scheduler();
    print_mmap_time();
}

uint64_t SKGgenerator::workload_size_calc_src_dim(int log_n, uint64_t src_vid_start, uint64_t src_vid_end){
    int log_postfix = get_lsb_loc(src_vid_end - src_vid_start);
    int bitcount = count_bits(src_vid_start >> log_postfix);
    double prob = pow(c+d, bitcount) * pow(a+b, log_n - log_postfix - bitcount);
    //assert(prob <= 1);
    uint64_t num_edge_in_workload =  static_cast<uint64_t>(n_edge * prob);
    return num_edge_in_workload;
}

uint64_t SKGgenerator::workload_size_calc_dst_dim(int log_n, uint64_t dst_vid_start, uint64_t dst_vid_end, uint64_t num_edge_in_src_dim){
    int log_postfix = get_lsb_loc(dst_vid_end - dst_vid_start);
    int bitcount = count_bits(dst_vid_start >> log_postfix);
    double prob = pow(b+d, bitcount) * pow(a+c, log_n - log_postfix - bitcount);
    uint64_t num_edge_in_workload =  static_cast<uint64_t>(num_edge_in_src_dim * prob);
    return num_edge_in_workload;
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


inline uint64_t flunctuate(uint64_t n_vertex, double prob){
    double stdev = sqrt(n_vertex * prob * (1 - prob));
    double mean = n_vertex * prob;

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(mean, stdev);

    uint64_t result = static_cast<uint64_t>(d(gen));
    return result;
}