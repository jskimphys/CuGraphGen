import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re

degree = np.zeros(1 << 22)

dir = 'graph'
for file in os.listdir(dir):
    if not re.match(r'part-\d+', file):
        continue
    num = np.loadtxt(dir + '/' + file, dtype=np.uint32, delimiter='\t')
    for i in range(num.shape[0]):
        degree[num[i, 0]] += 1
        degree[num[i, 1]] += 1

deg0_num = degree[degree == 0].shape[0]
print(f"deg0 vertex count : {deg0_num}")

x = degree[degree != 0]

logbins = np.logspace(np.log10(1),np.log10(1e6),60)
x.sort()
fig, ax = plt.subplots(figsize=(3,3))
ax.hist(x, bins=logbins )
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('degree')
ax.set_ylabel('count')
ax.set_ylim(1, 1e8)
ax.set_xlim(0.6, 1e6)
fig.tight_layout()
fig.savefig('trillgen_deg.png')