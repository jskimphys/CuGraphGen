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

## Disclaimer

**This project comes with no warranties or responsibilities.**

- **No Warranty:** CuGraphGen is provided as-is without any warranty. The software may have bugs, limitations, or inaccuracies, and its performance may vary.
  
- **No Maintenance Commitment:** The author may not actively maintain or update this project. Bug fixes, improvements, or support may not be provided.
  
- **Use at Your Own Risk:** Users are encouraged to carefully review the code and assess its suitability for their specific use case. The author is not responsible for any consequences resulting from the use or misuse of this software.
