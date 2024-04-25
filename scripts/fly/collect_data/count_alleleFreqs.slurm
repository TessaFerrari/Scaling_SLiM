#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --array=22,27-29,34-64,66-90 #90
#SBATCH --job-name=flyBench_sfs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/sfs/slurm-%A_%a.out

# Load necessary modules
module load gcc/9.2.0 samtools/18.0.4 htslib/1.10.2 bcftools/1.10.2 vcftools

# Set config directory
cd /scratch1/tferrari/SlimBenchmark/fly/scripts/data

# Select parameters from ArrayTaskID
config=fly_muts_config.txt
t=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)       # t=[5,10,20,Coal,Recap]
G=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)       # G=[1e5,1e6,1e7]
c=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)       # c=[1000,500,100]
n=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)       # n=[1,2]
h=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $6}' $config)       # h=[0.0,0.5]
rep=5

# Set dependent parameters
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/out/burn${t}_scale${c}_gensize${G}

mkdir sfs

# Collect mutation counts for each sim replicate
for pop in african; do

out=sfs/flyBench_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

#Print sfs header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tTotNumMuts" > ${out}
for i in {1..50}; do printf "\t$i"; done >> ${out}
printf "\n" >> ${out}

for (( r=1; r <= $z; r++ )); do

        file_bn=flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}.${pop}

        # Subset vcf to 50 random individuals and calc allele freqs 		| isolate allele freqs                      | print minor allele freqs                 | remove 0s
        freqs=`vcftools --gzvcf ${file_bn}.vcf.gz --max-indv 50 --freq --stdout | awk '{print substr($5,3), substr($6, 3)}' | awk '{if($1>$2)print $2; else print $1}' | awk '$1!=0 && NR!=1'`

        # Initialize sfs array with 0s
        declare -a SFS=( $(for i in {1..50}; do echo 0; done) )

        # Fill sfs using allele freqs (array index i-1 corresponds to freq i/100)
        for mut in $freqs; do
                (( SFS[$(echo "$mut*100-1" | bc | awk -F'.' '{print $1}')]++ ))
        done

        # Total number of snps
        totM=`echo "$freqs" | wc -l`

        printf "${pop}\t$t\t$c\t$G\t$h\t$n\t$r\t$totM" >> ${out}
        printf '\t%s' "${SFS[@]}" >> ${out}
	printf "\n" >> ${out}

done
done

