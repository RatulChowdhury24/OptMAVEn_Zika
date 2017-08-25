import os
for exp in os.listdir("."):
	if exp.startswith("bench"):
		os.system("qsub " + exp + "/start.sh")
