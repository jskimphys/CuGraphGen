#include <iostream>
#include <algorithm>
#include <string>
#include <filesystem>

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

    int edgefactor = 16;
    it = std::find(args.begin(), args.end(), "-e");
    if (it != args.end()) {
        edgefactor = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No edgefactor provided, using default value of 16" << std::endl;
    }

    uint64_t workload_size_limit = 1LL << 28; // 256MB
    it = std::find(args.begin(), args.end(), "-w");
    if (it != args.end()) {
        workload_size_limit = std::stoi(*(it + 1));
    }
    else{
        std::cout << "No workload size limit provided, using default value of 256MB" << std::endl;
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
    generator.schedule();
    generator.print_schedule();
}