import numpy as np
import matplotlib.pyplot as plt


filepath = 'out/deg.txt'

# [[out_degree, in_degree], ...]
array = np.loadtxt(filepath, dtype=np.uint32, delimiter='\t')
#find if there is a node with 0 in degree and 0 out degree
degree = array[:, 0] + array[:, 1]
print(np.where(degree == 0))

# plot the degree distribution

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