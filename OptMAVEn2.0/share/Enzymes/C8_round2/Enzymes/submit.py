# This script submits a list of numbered PBS scripts 

import os

num = 0
for i in range(10):
	num +=1
	job = "job_C8_round2_" + str(num)
	os.system("qsub " + job)

