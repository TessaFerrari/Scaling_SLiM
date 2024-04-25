#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --array=1-500 #1-1800 total
#SBATCH --job-name=flyBench_stats
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/stats/slurm-%A_%a.out

module load gcc/8.3.0
SLIMDIR=/home1/tferrari/build

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts

# Select parameters from ArrayTaskID
config=run/fly_main_config.txt 							# Use the same config as main
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)	# t=[5,10,20,Coal,Recap]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)	# G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)	# c=[10,5,1]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)	# n=[1,2]
h=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)       # h=[0.0,0.5]
rep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)     # rep=[1,2,3,4,5] or {if G=1e5}rep=[1..50]


# Generate and run  SLiM burn-in scripts (using different scripts if neutral mutations were overlaid or not)
if [[ "$t" == "Coal" ]]
then
        ./run/make_slim_fly_stats_over.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}				# generate script
        ${SLIMDIR}/slim -t -m ./temp/flyBench_stats_coal_c${c}_${G}_h${h}_n${n}_rep${rep}.job		# run simulation
elif [[ "$t" == "Recap" ]]
then
        ./run/make_slim_fly_stats_over.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}				# generate script
        ${SLIMDIR}/slim -t -m ./temp/flyBench_stats_recap_c${c}_${G}_h${h}_n${n}_rep${rep}.job		# run simulation
else
        ./run/make_slim_fly_stats_sim.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}				# generate script
        ${SLIMDIR}/slim -t -m ./temp/flyBench_stats_${t}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.job		# run simulation
fi


