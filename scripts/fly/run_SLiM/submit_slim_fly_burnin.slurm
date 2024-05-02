#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=70G
#SBATCH --array=1-72
#SBATCH --job-name=flyBench_burnin
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/slurm_logs/fly/burnin/slurm-%A_%a.out

module load gcc/8.3.0
SLIMDIR=/home1/tferrari/build

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/fly/run_SLiM

# Select parameters from config file
config=fly_burnin_config.txt 
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)	# t=[5,10,20,Coal]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)	# G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)	# c=[1000,500,100]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)	# n=[1,2]

# Create output directory
mkdir ../../../out/fly/burn${t}_scale${c}_gensize${G}

# Generate and run  SLiM burn-in scripts
if [[ "$t" == "Coal" ]]
then
	./slim_scripts/make_slim_fly_burnin_coal.job.sh ${n} ${c} ${G}		# generate coal script
	${SLIMDIR}/slim -t -m ../temp/flyBench_burnin${n}_coal_c${c}_${G}.job	# run coal simulation
else
	./slim_scripts/make_slim_fly_burnin_sim.job.sh ${n} ${t} ${c} ${G}		# generate 5N/10N/20N script
	${SLIMDIR}/slim -t -m ../temp/flyBench_burnin${n}_${t}Ne_c${c}_${G}.job		# run 5N/10N/20N simulation
fi

