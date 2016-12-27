import os

from matplotlib import pyplot as pl
import numpy as np

energies = np.array(sorted([float(open(f).readlines()[1].split()[4]) for f in os.listdir(os.getcwd()) if f.startswith("namdenergy") and f.endswith("dat")]))
pl.plot(np.log(np.abs(energies)))
pl.show()
