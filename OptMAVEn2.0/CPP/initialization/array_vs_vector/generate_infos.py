import numpy as np
for i, n in {1: 10000, 2: 1000}.items():
	open("coords{}.dat".format(i), "w").write('\n'.join([' '.join(["{:<5f}".format(x) for x in 10*np.random.random(3)]) for j in xrange(n)]))
