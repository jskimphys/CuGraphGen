#include <string>
#include <iostream>
#include <cstring>
#include <thread>
#include <cstdio>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <algorithm>

#include "fast_writer.h"
#include "constants.h"

using namespace std;

void mmap_write_DMA(string dst_filename, char* src_DMA_region, uint64_t size){
    int fd = open(dst_filename.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
    if (fd == -1) {
        perror("Error opening file for writing");
        exit(1);
    }
    if(ftruncate(fd, size) == -1){
        perror("Error truncating file");
        exit(1);
    }

    char* map = (char*)mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (map == MAP_FAILED) {
        close(fd);
        perror("Error mmapping the file");
        exit(1);
    }

    close(fd);

    memcpy(map, src_DMA_region, size);
    if (munmap(map, size) == -1) {
        perror("Error un-mmapping the file");
        exit(1);
    }
}

void fast_writer::write_async(string filepath,  char* mem_src, uint64_t size){
    while(writing){
        usleep(10);
    }
    auto job = [=](){
        writing = true;
        FILE *fp = fopen(filepath.c_str(), "w");
        fwrite(mem_src, sizeof(char), size, fp);
        // mmap_write_DMA(filepath, mem_src, size);
        writing = false;
    };
    thread t(job);
    usleep(10);
    t.detach();
}

void fast_writer::wait(){
    while(writing){
        usleep(10);
    }
}
