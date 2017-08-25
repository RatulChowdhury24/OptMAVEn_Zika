#PBS -l pmem=4gb
#PBS -l walltime=23:00:00
#PBS -j oe
#PBS -o /dev/null
#PBS -q lionxf

set -u
echo "Job started on $(hostname -s) at $(date)"
(time -p /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.sh) 2> /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.time.txt
cd /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/
python /gpfs/work/m/mfa5147/OptMAVEn2.0/modules/PERFORMANCE.py /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.time.txt MILP_for_positions_1036_to_1036
rm /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.sh /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.time.txt
rm /gpfs/work/m/mfa5147/OptMAVEn2.0/experiments/rep1_3bky_01/positions_1036_to_1036_MILP.timer.sh
echo "Job ended at $(date)"
