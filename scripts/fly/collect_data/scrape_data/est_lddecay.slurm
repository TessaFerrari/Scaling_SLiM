#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --array=1-90 #90
#SBATCH --job-name=flyBench_ld
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/ld/slurm-%A_%a.out

# Load necessary modules
module load vcftools

# Set config directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts/data

# Set dist, the max distance between snps for binning
dist=1000

# Set nbins, the number of bins to create
nbins=10

# Select parameters from ArrayTaskID
config=fly_muts_config.txt
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)       # t=[5,10,20,Coal,Recap]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)       # G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)       # c=[10,5,1]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)       # n=[1,2]
h=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)       # h=[0.0,0.5]
rep=5

# Set dependent parameters
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/out/burn${t}_scale${c}_gensize${G}

mkdir ld

# Collect mutation counts for each sim replicate
for pop in african; do

out=ld/dist${dist}_${nbins}bins_flyBench_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

#Print ld header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum" > ${out}
for (( i=1; i <= $nbins; i++ )); do printf "\tb"$i"_nSNPs\tb"$i"_sumR2\tb"$i"_sumD\tb"$i"_Dp"; done >> ${out}
printf "\n" >> ${out}

	# Loop through each simulation replicate
	for (( r=1; r <= $z; r++ )); do

        	file_bn=flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}.${pop}

		# output vcftools ld stats
		ld_stats=`vcftools --gzvcf ${file_bn}.vcf.gz --hap-r2 --ld-window-bp ${dist} --maf 0.1 --stdout | awk 'BEGIN{OFS="\t"} {print $3-$2, $5, $6, $7}'`

		# print sim ID stats
		printf "${pop}\t$t\t$c\t$G\t$h\t$n\t$r" >> ${out}
	
		# print LD stats for each bin
		for (( bin=0; bin < $nbins; bin++ )) ; do

			printf "$ld_stats" | awk -v b=$bin -v x=$nbins -v y=$dist 'BEGIN{sumR2=0; sumD=0; sumDp=0; i=0} $1>=(b*(y/x)) && $1<((b+1)*(y/x)){sumR2+=$2; sumD+=$3; sumDp+=$4; i++} END{OFS="\t";printf "\t"i"\t"sumR2"\t"sumD"\t"sumDp}'

		done >> ${out}
        
		# newline 
		printf "\n" >> ${out}

	done
done

