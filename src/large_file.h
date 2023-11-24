
#include <filesystem>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <iostream>
#include <chrono>

using namespace std;

uint64_t get_avail_storage(std::string& dir){
    return std::filesystem::space(dir).available;
}

chrono::milliseconds open_time;
chrono::milliseconds mmap_time;
chrono::milliseconds munmap_time;

char* get_mmap_memory(string filepath,  uint64_t size){
    auto start = std::chrono::high_resolution_clock::now();
    
    
    //mmap to filepath
    start = std::chrono::high_resolution_clock::now();
    int fd = open(filepath.c_str(), O_RDWR | O_CREAT, 0666);
    if (fd == -1){
        std::cout << "open file failed" << std::endl;
        exit(0);
    }
    if (ftruncate(fd, size) == -1){//extend file to size
        std::cout << "ftruncate failed" << std::endl;
        exit(0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    open_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    start = std::chrono::high_resolution_clock::now();
    char* mem = (char*)mmap(NULL, size, PROT_READ | PROT_WRITE , MAP_SHARED, fd, 0);
    if (mem == MAP_FAILED){
        perror("mmap failed");
        exit(0);
    }
    close(fd);
    end = std::chrono::high_resolution_clock::now();
    mmap_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return mem;
}

void free_mmap_memory(char* mem, uint64_t size){
    //shrunk memorysize
    auto start = std::chrono::high_resolution_clock::now();
    if (munmap(mem, size) == -1){
        std::cout << "munmap failed" << std::endl;
        exit(0);
    }
    auto end = std::chrono::high_resolution_clock::now();
    munmap_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
}

void print_mmap_time(){
    cout << "open time: " << open_time.count() << "ms" << endl;
    cout << "mmap time: " << mmap_time.count() << "ms" << endl;
    cout << "munmap time: " << munmap_time.count() << "ms" << endl;
}

