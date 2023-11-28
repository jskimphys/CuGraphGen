
CC=nvcc
sm_arch=sm_75
sm_flags=--generate-code arch=compute_75,code=sm_75
opt_flags=-O3
lib_flags=-lcuda -lcurand
flags=$(sm_flags) $(opt_flags) $(lib_flags)

all: main test

main: src/*.cu src/*.cpp
	$(CC) $(flags) -o $@ $^

test: util/*.cpp
	g++ -O2 -o $@ $^


clean:
	rm -f main