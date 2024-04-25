#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --array=1 #90
#SBATCH --job-name=flyBench_sel
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tferrari@usc.edu
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --output=/scratch1/tferrari/SlimBenchmark/fly/logs/slurm/mutStats/slurm-%A_%a.out

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

# Set working directory
cd /scratch1/tferrari/SlimBenchmark/fly/out/burn${t}_scale${c}_gensize${G}

mkdir mutStats

# Collect mutation counts for each sim replicate
for pop in african; do

out=mutStats/flyBench_start_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

#Print header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tTot.MeanS\tDel.MeanS\tTot.MeanTO\tDel.MeanTO\tNeu.MeanTO\tCurrentTick\tStartingTick\n" > ${out}
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tStartingTick\n" > ${out}

for (( r=1; r <= $z; r++ )); do

        file_bn=flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}

        #totS=`zcat ${file_bn}.${pop}.vcf.gz | awk 'NR>14{print $8}' | awk -F';' '{print substr($2, 3)}' | awk '{ sumS += $1; n++ } END { if (n > 0) print sumS/n; }'`
	#delS=`zcat ${file_bn}.${pop}.vcf.gz | awk 'NR>14{print $8}' | awk -F';' 'substr($2, 3)<0{print substr($2, 3)}' | awk '{ sumS += $1; n++ } END { if (n > 0) print sumS/n; }'`
	#totTO=`zcat ${file_bn}.${pop}.vcf.gz | awk 'NR>14{print $8}' | awk -F';' '{print substr($5, 4)}' | awk '{ sumT += $1; n++ } END { if (n > 0) print sumT/n; }'`
	#delTO=`zcat ${file_bn}.${pop}.vcf.gz | awk 'NR>14{print $8}' | awk -F';' 'substr($2, 3)<0{print substr($5, 4)}' | awk '{ sumT += $1; n++ } END { if (n > 0) print sumT/n; }'`
	#neuTO=`zcat ${file_bn}.${pop}.vcf.gz | awk 'NR>14{print $8}' | awk -F';' 'substr($2, 3)==0{print substr($5, 4)}' | awk '{ sumT += $1; n++ } END { if (n > 0) print sumT/n; }'`
	#CT=`awk 'NR==37{print $1}' ../../scripts/temp/${file_bn}.job`
	ST=`awk 'NR==23{print $1}' ../../scripts/temp/${file_bn}.job`

        #printf "$pop\t$t\t$c\t$G\t$h\t$n\t$r\t$totS\t$delS\t$totTO\t$delTO\t$neuTO\t$CT\t$ST\n" >> ${out}
	printf "$pop\t$t\t$c\t$G\t$h\t$n\t$r\t$ST\n" >> ${out}
done
done

