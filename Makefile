
CC=nvcc
sm_arch=sm_75
sm_flags=--generate-code arch=compute_75,code=sm_75
opt_flags=-O3
lib_flags=-lcuda -lcurand
flags=$(sm_flags) $(opt_flags) $(lib_flags)

all: main

main: src/*.cu src/*.cpp
	$(CC) $(flags) -o main src/*.cu src/*.cpp 

test:


clean:
	rm -f main