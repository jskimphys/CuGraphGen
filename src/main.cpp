#include <iostream>
#include <algorithm>
#include <string>
#include <filesystem>
#include <chrono>

#include "SKGgenerator.h"


int main(int argc, char** argv) {
    //parse string
    std::vector<std::string> args(argv, argv + argc);

    int scale = 20;
    auto it = std::find(args.begin(), args.end(), "-s");
    if (it != args.end()) {
        scale = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No scale provided, using default value of 20" << std::endl;
    }

    int edgefactor = 10;
    it = std::find(args.begin(), args.end(), "-e");
    if (it != args.end()) {
        edgefactor = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No edgefactor provided, using default value of 10" << std::endl;
    }

    uint64_t workload_size_limit = 1ULL << 26; // 64MB
    it = std::find(args.begin(), args.end(), "-w");
    if (it != args.end()) {
        workload_size_limit = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No workload size limit provided, using default value of 64MB" << std::endl;
    }

    int seed = 0;
    it = std::find(args.begin(), args.end(), "--seed");
    if (it != args.end()) {
        seed = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No seed provided, using default value of 0" << std::endl;
    }

    /*
    //make directory to save graphs
    std::string dir = "out";
    it = std::find(args.begin(), args.end(), "-d");
    if (it != args.end()) {
        dir = *(it + 1);
    }
    else{
        std::cout << "No directory provided, using default value of out" << std::endl;
    }

    //generate directory
    //directory already exists
    if (std::filesystem::exists(dir)){
        std::cout << "Directory already exists, exiting" << std::endl;
        return 0;
    }
    else{
        std::filesystem::create_directory(dir);
    }
    */

    //generate graphs
    SKGgenerator generator(0.57, 0.19, 0.19, 0.05, scale, scale + edgefactor, seed, "out", workload_size_limit);
    auto start = std::chrono::high_resolution_clock::now();
    generator.schedule();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time to generate schedule: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    generator.print_workload_summary();
    //generator.print_workload();


    generator.generate();
}