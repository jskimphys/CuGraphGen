#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>

using namespace std;

int main(){
    ifstream fs("out/file_info.json");
    int n_file = 0;
    uint32_t n_vertex = 0;
    uint64_t n_edge = 0;
    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;

    //{
// 	"n_file" : 51,
// 	"n_vertex" : 134217728,
// 	"n_edge" : 2147483648,
// 	"a" : 0.57,
// 	"b" : 0.19,
// 	"c" : 0.19,
// 	"d" : 0.05
// }
    string str;
    while(fs >> str){
        if(str.find("n_file") != string::npos){
            fs >> str;
            fs >> str;
            n_file = stoi(str.substr(0, str.size()-1));
        }
        else if(str.find("n_vertex") != string::npos){
            fs >> str;
            fs >> str;
            n_vertex = stoi(str.substr(0, str.size()-1));
        }
        else if(str.find("n_edge") != string::npos){
            fs >> str;
            fs >> str;
            n_edge = stoull(str.substr(0, str.size()-1));
        }
        else if(str.find("a") != string::npos){
            fs >> str;
            fs >> str;
            a = stod(str.substr(0, str.size()-1));
        }
        else if(str.find("b") != string::npos){
            fs >> str;
            fs >> str;
            b = stod(str.substr(0, str.size()-1));
        }
        else if(str.find("c") != string::npos){
            fs >> str;
            fs >> str;
            c = stod(str.substr(0, str.size()-1));
        }
        else if(str.find("d") != string::npos){
            fs >> str;
            fs >> str;
            d = stod(str.substr(0, str.size()-1));
        }
    }
    cout << n_file << "\t" << n_vertex << "\t" << n_edge << "\t" << a << "\t" << b << "\t" << c << "\t" << d << endl;

    vector<uint32_t> out_deg(n_vertex, 0);
    vector<uint32_t> in_deg(n_vertex, 0);

    for(int i=0; i<n_file; i++){
        string file_idx = to_string(i);
        while(file_idx.size() < 3){
            file_idx = "0" + file_idx;
        }
        string path = "out/edgelist8B_" + file_idx + ".part";
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
        if(mem == MAP_FAILED){
            std::cout << "mmap failed" << std::endl;
            exit(0);
        }

        uint64_t* mem64 = (uint64_t*)mem;
        for(int i=0; i<sb.st_size/16; i++){
            out_deg[mem64[2*i]]++;
            in_deg[mem64[2*i+1]]++;
        }
        cout << "finish file " << i << endl;
        munmap(mem, sb.st_size);
        close(fd);
    }

    string deg_file_path = "out/deg.txt";
    FILE *deg_file = fopen(deg_file_path.c_str(), "w");
    for(int i=0; i<n_vertex; i++){
        fprintf(deg_file, "%u\t%u\n", out_deg[i], in_deg[i]);
    }
}
