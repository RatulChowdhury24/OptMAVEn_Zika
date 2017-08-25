#PBS -l pmem=4gb
#PBS -l walltime=23:59:59
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2b1h_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2b1h_01
cd rep1_2b1h_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2b1h_antigen.pdb 2b1h_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2b1h_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2iff_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2iff_01
cd rep1_2iff_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2iff_antigen.pdb 2iff_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2iff_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2a6i_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2a6i_01
cd rep1_2a6i_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2a6i_antigen.pdb 2a6i_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2a6i_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2zpk_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_2zpk_01
cd rep1_2zpk_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2zpk_antigen.pdb 2zpk_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/2zpk_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1kc5_01
cd rep1_1kc5_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kc5_antigen.pdb 1kc5_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kc5_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1kc5_inputs2.txt
