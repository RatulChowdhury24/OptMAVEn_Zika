import sys

import numpy as np

with open("parameters.txt") as p:
	line = p.readline()
	while not line.startswith("Epitope"):
		line = p.readline()
	epitope = set(map(int, line.split()[1:]))

allcom = np.mean(np.array([map(float, [line[31: 38], line[39: 46], line[47: 54]]) for line in open(sys.argv[1]).readlines() if line.startswith("ATOM")]), axis=0)

epicom = np.mean(np.array([map(float, [line[31: 38], line[39: 46], line[47: 54]]) for line in open(sys.argv[1]).readlines() if line.startswith("ATOM") and int(line[23: 27]) in epitope]), axis=0)

print "All COM:", allcom
print "Epi COM:", epicom
