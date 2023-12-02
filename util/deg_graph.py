import numpy as np
import matplotlib.pyplot as plt


filepath = 'out/deg.txt'

# [[out_degree, in_degree], ...]
array = np.loadtxt(filepath, dtype=np.uint32, delimiter='\t')
#find if there is a node with 0 in degree and 0 out degree
degree = array[:, 0] + array[:, 1]

# plot the degree distribution

deg0_num = degree[degree == 0].shape[0]
print(f"deg0 vertex count : {deg0_num}")

x = degree[degree != 0]
# print(bins)
logbins = np.logspace(np.log10(1),np.log10(1e6),60)
x.sort()
print(x[x<1e4][-1])
print(x[x>1e4][0])
#setting figure for the paper
fig, ax = plt.subplots(figsize=(3,3))
ax.hist(x, bins=logbins )
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('degree')
ax.set_ylabel('count')
ax.set_ylim(1, 1e8)
ax.set_xlim(0.6, 1e6)
fig.tight_layout()
fig.savefig('cuGen_deg.png')