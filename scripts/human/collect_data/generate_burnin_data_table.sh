
cd /scratch1/tferrari/SlimBenchmark/human/out

printf "BurnInType\tScalingFactor\tGenomeSize\tRepNum\tMemMB\tTimeSecs\n" > ../data/burnin_compStat_table.txt
printf "ScalingFactor\tGenomeSize\tRepNum\tCoalTime\n" > ../data/burnin_coalescenceTimes_table.txt

for t in 5 10 20 Coal; do
for c in 10 5 1; do
for G in 1e5 1e6 1e7; do
for n in 1 2; do

	if [[ "$t" != "Coal" ]]; then x=$t"Ne"; else x=coal; fi
	file=burn${t}_scale${c}_gensize${G}/humanBench_burnin${n}_${x}_c${c}_${G}.txt

	# Computational Efficiency Stats
	mem=`awk -F", |MB)" '/Peak memory usage:/ {print $2}' $file`
	time=`awk -F": " '/CPU time used:/ {print $2}' $file`

	printf "$t\t$c\t$G\t$n\t$mem\t$time\n" >> ../data/burnin_compStat_table.txt

	# Coalescence Stats
	if [[ "$t" == "Coal" ]]; then
		coal=`awk -F"[][()]" '/COALESCE/ {print $2}' $file`
		printf "$c\t$G\t$n\t$coal\n" >> ../data/burnin_coalescenceTimes_table.txt
	fi


done
done
done
done
