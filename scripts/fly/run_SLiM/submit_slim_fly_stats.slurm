#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --array=1-1800
#SBATCH --job-name=flyBench_stats
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/slurm_logs/fly/stats/slurm-%A_%a.out

module load gcc/8.3.0
SLIMDIR=/home1/tferrari/build

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/fly/run_SLiM

# Select parameters from config file (same as main)
config=fly_main_config.txt 
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)	# t=[5,10,20,Coal,Recap]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)	# G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)	# c=[1000,500,100]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)	# n=[1,2]
h=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)	# h=[0.5]
rep=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $config)	# rep=[1,2,3,4,5] or {if G=1e5}rep=[1..50]


# Generate and run  SLiM burn-in scripts
if [[ "$t" == "Coal" ]]
then
	./slim_scripts/make_slim_fly_stats_over.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}		# generate coal script
	${SLIMDIR}/slim -t -m ../temp/flyBench_stats_coal_c${c}_${G}_h${h}_n${n}_rep${rep}.job	# run coal simulation
elif [[ "$t" == "Recap" ]]
then

	./slim_scripts/make_slim_fly_stats_over.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}		# generate recap script
	${SLIMDIR}/slim -t -m ../temp/flyBench_stats_recap_c${c}_${G}_h${h}_n${n}_rep${rep}.job	# run recap simulation
else
	./slim_scripts/make_slim_fly_stats_sim.job.sh ${n} ${t} ${c} ${G} ${h} ${rep}			# generate 5N/10N/20N script
	${SLIMDIR}/slim -t -m ../temp/flyBench_stats_${t}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.job	# run 5N/10N/20N simulation
fi

