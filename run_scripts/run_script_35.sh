#PBS -l pmem=4gb
#PBS -l walltime=23:59:59
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3iu3_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_3iu3_01
cd rep1_3iu3_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3iu3_antigen.pdb 3iu3_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3iu3_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kcr_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1kcr_01
cd rep1_1kcr_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kcr_antigen.pdb 1kcr_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kcr_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3qg6_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_3qg6_01
cd rep1_3qg6_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3qg6_antigen.pdb 3qg6_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3qg6_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2r29_01
cd rep1_2r29_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2r29_antigen.pdb 2r29_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2r29_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2r29_inputs2.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1ggi_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1ggi_01
cd rep1_1ggi_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1ggi_antigen.pdb 1ggi_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1ggi_inputs1.txt
