# What is this program?
This is an simple proof of concept program that runs R-MAT graph generation algorithm in CUDA GPU.

# how to run?
First, in the Makefile, change the sm_arch that matches your GPU architecture.
You also need cuda toolkit installed in your system.

Then, run:
```
make
```
will simply compile the program.

the options are
-s: scale (|V| = 2^s)
-e: edgefactor (|E| = e * |V|)
--seed: seed for random number generator
-d: directory to save the result
--w:workload size limit(ignoring is okay unless the prgram fails)

ex) ./main -s 25

Each options will be set to default if not specified. and print the defualt values.

# result format
result will be saved in the directory specified by -d option.
result file is is a binary file edgelist, where edge (u, v) is represented by two 32-bit unsigned integers u and v.
If 32bit is not enough to represent the vertex id, you can change "src/constants.h" file.
