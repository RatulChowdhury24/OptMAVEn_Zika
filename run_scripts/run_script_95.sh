#PBS -l pmem=4gb
#PBS -l walltime=23:59:59
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1xgy_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1xgy_01
cd rep1_1xgy_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1xgy_antigen.pdb 1xgy_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1xgy_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1obe_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1obe_01
cd rep1_1obe_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1obe_antigen.pdb 1obe_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1obe_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1pz5_01
cd rep1_1pz5_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1pz5_antigen.pdb 1pz5_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1pz5_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1pz5_inputs2.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1acy_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1acy_01
cd rep1_1acy_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1acy_antigen.pdb 1acy_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1acy_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1bj1_01
cd rep1_1bj1_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1bj1_antigen.pdb 1bj1_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1bj1_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1bj1_inputs2.txt
