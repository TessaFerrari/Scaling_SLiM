#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --array=1-64,66-90 #90
#SBATCH --job-name=flyBench_muts
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/muts/slurm-%A_%a.out

# Load necessary modules
module load gcc/9.2.0 samtools/18.0.4 htslib/1.10.2 bcftools/1.10.2

# Set config directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts/data

# Select parameters from ArrayTaskID
config=fly_muts_config.txt
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)	# t=[5,10,20,Coal,Recap]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)	# G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)	# c=[10,5,1]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)	# n=[1,2]
h=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)       # h=[0.0,0.5]
rep=5

# Set dependent parameters
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/out/burn${t}_scale${c}_gensize${G}

mkdir mut_counts

# Collect mutation counts for each sim replicate
for (( r=1; r <= $z; r++ )); do

	# Set file basename
	file_bn=flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}

	# Collect mutation counts for each pop
	for pop in african; do

		# Format vcf file
		bgzip ${file_bn}.${pop}.vcf
		tabix -p vcf ${file_bn}.${pop}.vcf.gz

		# Collect deleterious and neutral counts
		del=`gunzip -c ${file_bn}.${pop}.vcf.gz | grep -o -i S=- | wc -l`
		neu=`gunzip -c ${file_bn}.${pop}.vcf.gz | grep -o -i S=0 | wc -l`

		# Output counts to file
		printf "$pop deleterious: $del\n" >> mut_counts/${file_bn}.mutcounts.txt
		printf "$pop neutral: $neu\n" >> mut_counts/${file_bn}.mutcounts.txt
	done
done

