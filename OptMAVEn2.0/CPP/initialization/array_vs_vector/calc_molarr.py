import os
import sys

cpp = ["""#include <fstream>
#include "molarr.cpp"
using namespace std;

int main()
{
	MolArr mol1, mol2;"""]

for i in range(1, 3):
	cpp.append("	float coords[][3] = {{{}}};".format(", ".join(["{{{}}}".format(", ".join(line.split())) for m, line in enumerate(open(sys.argv[i]))])))
	cpp.append("	mol{}.size = {};".format(i, m + 1))
	cpp.append("""	for (int i = 0; i < {}; i++)
			{
				for 
				coords
			}
		""".format(m)

	

cpp.append("""
	return 0;
}
""")

code = "calc_molarr"
src = code + ".cpp"
out = code + ".out"
open(src, "w").write("\n".join(cpp))
#os.system("g++ {} -o {}".format(src, out))
