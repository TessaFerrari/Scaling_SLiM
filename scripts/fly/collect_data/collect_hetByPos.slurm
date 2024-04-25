#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --array=1 #90
#SBATCH --job-name=flyBench_het
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/het/slurm-%A_%a.out

# Load necessary modules
module load gcc/9.2.0 samtools/18.0.4 htslib/1.10.2 bcftools/1.10.2 vcftools

# Set config directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts/data

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

# Set window size
nWins=20
if [[ "$G" == "1e5" ]]; then len=100000; elif [[ "$G" == "1e6" ]]; then len=1000000; else len=10000000; fi
winSize=$(($len/$nWins))

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/out/burn${t}_scale${c}_gensize${G}

mkdir het

# Collect mutation counts for each sim replicate
for pop in african; do

out=het/flyBench_start_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

#Print header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum" > ${out}
for (( w=1; w <= $nWins; w++ )); do printf "\tneu_$i"; done >> ${out}
for (( w=1; w <= $nWins; w++ )); do printf "\tdel_$i"; done >> ${out}
for (( w=1; w <= $nWins; w++ )); do printf "\thetsum_$i"; done >> ${out}
printf "\n" >> ${out}

for (( r=1; r <= $z; r++ )); do

        file_bn=flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}

	zcat ${file_bn}.${pop}.vcf.gz | awk -v OFS="\t" 'NR>14{ for (i=10;i<=NF;i++) if ($i=="0|1" || $i=="1|0") c++; split($8,s,";"); print $2, substr(s[2], 3), c/(NF-9); c=0 }' > het/tempData.txt

	# Initialize array with size $nWins
        declare -a NEU=( $(for (( w=1; w <= $nWins; w++ )); do echo 0; done) ) # counts neutral mutations for each window
	declare -a DEL=( $(for (( w=1; w <= $nWins; w++ )); do echo 0; done) ) # counts deleterious mutations for each window
	declare -a HETSUM=( $(for (( w=1; w <= $nWins; w++ )); do echo 0; done) ) # sums heterozygosity for each window

        # Fill in array using mutations
	while read line; do
		pos=`echo "$line" | awk '{print $1}'`
		S=`echo "$line" | awk '{print $2}'`
		het=`echo "$line" | awk '{print $3}'`

		windex=$(( $pos/$winSize ))

                if [[ $S==0 ]]; then (( NEU[$windex]++ ))
        	else (( DEL[$windex]++ )); fi
		currHetSum=${HETSUM[$windex]}
		HETSUM[$windex]=$(echo "scale=2; $currHetSum + $het" | bc)
	done < het/tempData.txt

	printf "$pop\t$t\t$c\t$G\t$h\t$n\t$r" >> ${out}
	printf '\t%s' "${NEU[@]}" >> ${out}
	printf '\t%s' "${DEL[@]}" >> ${out}
	printf '\t%s' "${HETSUM[@]}" >> ${out}
	printf "\n" >> ${out}

done
done

