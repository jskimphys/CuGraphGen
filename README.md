
cuda_sm=75
cuda_arch=compute_${cuda_sm},code=sm_${cuda_sm}


all: main test

main: main
    nvcc -arch=${cuda_arch}  src/*.cu src/*.cpp -lcuda -lcurand -O3 -o main

test: test

clean:
    rm -f main
