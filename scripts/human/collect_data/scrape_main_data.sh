
cd /scratch1/tferrari/SlimBenchmark/human

set -o noclobber # prevent overwriting if there are duplicate slurm logs 

for file in logs/slurm/main/*
do
	# Collect info for file naming
	n=`awk '/Burn-in replicate number:/ {print $6}' $file`
	t=`awk '/Burn-in type:/ {print $5}' $file`
	c=`awk '/Scaling factor:/ {print $5}' $file`
	G=`awk '/Genome size:/ {print $5}' $file`
	h=`awk '/Dominance coefficient:/ {print $5}' $file`
	rep=`awk '/Simulation replicate number:/ {print $6}' $file`

	# Add Ne suffix to numeric t values
	if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

	# Move stats to data file
	awk '$2 ~ /\*\*\*\*\*\*\*\*\*\*/ {print $0}' $file | cut -d\  -f3- > out/burn${t}_scale${c}_gensize${G}/humanBench_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.txt # print data to file

done

