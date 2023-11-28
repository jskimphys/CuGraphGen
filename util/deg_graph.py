import numpy as np
import matplotlib.pyplot as plt


filepath = 'out/deg.txt'

# [[out_degree, in_degree], ...]
array = np.loadtxt(filepath, dtype=np.uint32, delimiter='\t')
#find if there is a node with 0 in degree and 0 out degree
degree = array[:, 0] + array[:, 1]
print(np.where(degree == 0))

exit(0)
plt.hist(np.log(degree), bins=100)
plt.xscale('log')
plt.xlabel('out degree')
plt.ylabel('coun')
plt.show()