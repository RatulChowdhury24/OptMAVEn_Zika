#PBS -l pmem=4gb
#PBS -l walltime=23:59:59
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf
set -u

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzi_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1tzi_01
cd rep1_1tzi_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzi_antigen.pdb 1tzi_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzi_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1mvu_01
cd rep1_1mvu_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1mvu_antigen.pdb 1mvu_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1mvu_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1mvu_inputs2.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1i9r_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1i9r_01
cd rep1_1i9r_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1i9r_antigen.pdb 1i9r_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1i9r_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzg_inputs2.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_1tzg_01
cd rep1_1tzg_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzg_antigen.pdb 1tzg_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/1tzg_inputs1.txt



cd /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0
mkdir rep1_3d85_01
cd rep1_3d85_01
cp /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3d85_antigen.pdb 3d85_antigen.pdb
python /gpfs/scratch/mfa5147/BENCHMARKING_1.0_vs_2.0/OptMAVEn1.0/scripts/start_optmaven.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3d85_inputs1.txt

cd /gpfs/home/mfa5147/work/OptMAVEn2.0
python Start_Experiment.py < /gpfs/home/mfa5147/work/OptMAVEn2.0/input_files/3d85_inputs2.txt
