import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re

degree = np.zeros(1 << 22)

dir = 'graph_trillion'
for file in os.listdir(dir):
    if not re.match(r'part-\d+', file):
        continue
    num = np.loadtxt(dir + '/' + file, dtype=np.uint32, delimiter='\t')
    for line in open(dir + '/' + file):
        src, dst = line.split()
        degree[int(src)] += 1
        degree[int(dst)] += 1


x = degree[degree != 0]
hist, bins = np.histogram(x, bins=100)
logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
x.sort()
print(x[x<1e4][-1])
print(x[x>1e4][0])
plt.hist(x, bins=logbins )
plt.xscale('log')
plt.yscale('log')
plt.ylabel('out degree')
plt.xlabel('count')
plt.ylim(1, 1e8)
plt.show()