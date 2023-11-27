#include <iostream>
#include <algorithm>
#include <string>
#include <filesystem>
#include <chrono>

#include "SKGgenerator.cuh"

const int DEFAULT_SCALE = 20;
const int DEFAULT_EDGERATIO = 16;
const int DEFAULT_WORKLOAD_SIZE_LIMIT = 1ULL << 26;
const int DEFAULT_SEED = 0;
const std::string DEFAULT_DIR = "out";

const double a = 0.57;
const double b = 0.19;
const double c = 0.19;
const double d = 0.05;


int main(int argc, char** argv) {
    //parse string
    std::vector<std::string> args(argv, argv + argc);

    int scale = DEFAULT_SCALE;
    auto it = std::find(args.begin(), args.end(), "-s");
    if (it != args.end()) {
        scale = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No scale provided, using default value of " << DEFAULT_SCALE << std::endl;
    }

    uint64_t edgeratio = DEFAULT_EDGERATIO;
    it = std::find(args.begin(), args.end(), "-e");
    if (it != args.end()) {
        edgeratio = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No edge ratio provided, using default value of " << DEFAULT_EDGERATIO << std::endl;
    }

    uint64_t workload_byte_limit = DEFAULT_WORKLOAD_SIZE_LIMIT;
    it = std::find(args.begin(), args.end(), "-w");
    if (it != args.end()) {
        workload_byte_limit = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No workload size limit provided, using default value of " << DEFAULT_WORKLOAD_SIZE_LIMIT << std::endl;
    }

    int seed = DEFAULT_SEED;
    it = std::find(args.begin(), args.end(), "--seed");
    if (it != args.end()) {
        seed = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No seed provided, using default value of " << DEFAULT_SEED << std::endl;
    }


    std::string dir = DEFAULT_DIR;
    it = std::find(args.begin(), args.end(), "-d");
    if (it != args.end()) {
        dir = *(it + 1);
    }
    else{
        std::cout << "No directory provided, using default value of " << DEFAULT_DIR << std::endl;
    }

    //generate graphs
    SKGgenerator generator(a, b, c, d, scale, edgeratio, seed, dir, workload_byte_limit);
    auto start = std::chrono::high_resolution_clock::now();
    generator.divide_workloads();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time to generate divide_workloads: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    generator.print_workload_summary();

    start = std::chrono::high_resolution_clock::now();
    generator.generate();
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Time to generate graphs: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
    std::cout << "successfully generated graphs" << std::endl;
}