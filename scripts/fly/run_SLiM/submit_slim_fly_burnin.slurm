#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --mem=70G
#SBATCH --array=11-18,29-36,47-54,66-72 #72
#SBATCH --job-name=flyBench
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/burnin/slurm-%A_%a.out

module load gcc/8.3.0
SLIMDIR=/home1/tferrari/build

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts

# Select parameters from ArrayTaskID
config=run/fly_burnin_config.txt 
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)	# t=[5,10,20,Coal]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)	# G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)	# c=[1000,100,10]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)	# n=[1,2]

# Create output directory
mkdir ../out/burn${t}_scale${c}_gensize${G}

# Generate and run  SLiM burn-in scripts (using different scripts if burn-in is running until coalesced vs 5/10/20Ne generations)
if [[ "$t" == "Coal" ]]
then
	./run/make_slim_fly_burnin_coal.job.sh ${n} ${c} ${G}			# generate script
	${SLIMDIR}/slim -t -m ./temp/flyBench_burnin${n}_coal_c${c}_${G}.job	# run simulation
else
	./run/make_slim_fly_burnin_sim.job.sh ${n} ${t} ${c} ${G}		# generate script
	${SLIMDIR}/slim -t -m ./temp/flyBench_burnin${n}_${t}Ne_c${c}_${G}.job	# run simulation
fi

