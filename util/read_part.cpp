#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <string>

using namespace std;

int main(){
    string path = "../src/out/edgelist8B_0.part";
    int fd = open(path.c_str(), O_RDONLY);
    if (fd == -1){
        std::cout << "open file failed" << std::endl;
        exit(0);
    }

    struct stat sb;
    if (fstat(fd, &sb) == -1){
        std::cout << "fstat failed" << std::endl;
        exit(0);
    }

    uint64_t* mem = (uint64_t*)mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

    for(int i=0; i<20; i++){
        cout << mem[2*i] << " " << mem[2*i+1] << endl;
    }
}