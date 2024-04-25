
cd /scratch1/tferrari/SlimBenchmark/fly

set -o noclobber # prevent overwriting if there are duplicate slurm logs 

for file in logs/slurm/burnin/*
do
	# Collect info for file naming
	n=`awk '/Burn-in replicate number:/ {print $6}' $file`
	t=`awk '/Burn-in type:/ {print $5}' $file`
	c=`awk '/Scaling factor:/ {print $5}' $file`
	G=`awk '/Genome size:/ {print $5}' $file`

	# Add Ne suffix to numeric t values
	if [[ "$t" != "Coal" ]]; then x=$t"Ne"; else x=coal; fi

	# Move stats to data file
	awk '$2 ~ /\*\*\*\*\*\*\*\*\*\*/ {print $0}' $file | cut -d\  -f3- > out/burn${t}_scale${c}_gensize${G}/flyBench_burnin${n}_${x}_c${c}_${G}.txt # print data to file

	# Append coalescence info to file
	awk '/COALESCE/ {print $0}' $file >> out/burn${t}_scale${c}_gensize${G}/flyBench_burnin${n}_${x}_c${c}_${G}.txt

done

