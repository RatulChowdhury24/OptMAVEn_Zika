#PBS -l pmem=4gb
#PBS -l walltime=23:59:59
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2ck0_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2ck0_01
cd rep1_2ck0_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2ck0_antigen.pdb 2ck0_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2ck0_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1jhl_01
cd rep1_1jhl_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1jhl_antigen.pdb 1jhl_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1jhl_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1jhl_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2osl_01
cd rep1_2osl_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2osl_antigen.pdb 2osl_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2osl_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2osl_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_3eo1_01
cd rep1_3eo1_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3eo1_antigen.pdb 3eo1_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3eo1_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3eo1_inputs2.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3nfp_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_3nfp_01
cd rep1_3nfp_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3nfp_antigen.pdb 3nfp_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3nfp_inputs1.txt
