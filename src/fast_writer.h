#pragma once
#include <string>
#include <vector>

using namespace std;

class fast_writer{
public:
    fast_writer() :  writing(false){}
    void write_async(string filepath,  char* mem_src, uint64_t size);
    void wait();
private:
    volatile bool writing;
};