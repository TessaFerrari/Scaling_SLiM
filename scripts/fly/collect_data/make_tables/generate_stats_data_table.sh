# Set number of replicates per burn-in
rep=5 

cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/fly

printf "BurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tAHet\tATheta\n" > ../../data/fly/stats_popStat_table.txt
printf "BurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tMemMB\tTimeSecs\n" > ../../data/fly/stats_compStat_table.txt

for t in 5 10 20 Coal Recap; do
for c in 1000 500 100; do
for G in 1e5 1e6 1e7; do
for h in 0.5; do
for n in 1 2; do
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi
for (( r=1; r <= $z; r++ )); do

	if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi
	file=burn${t}_scale${c}_gensize${G}/flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}.stats.txt

	# Population stats
	a_het=`awk -F": " '/African heterozygousity:/ {print $2}' $file`
	a_theta=`awk -F": " '/African theta:/ {print $2}' $file`

	# Computational Efficiency Stats
	mem=`awk -F", |MB)" '/Peak memory usage:/ {print $2}' $file`
	time=`awk -F": " '/CPU time used:/ {print $2}' $file`

	printf "$t\t$c\t$G\t$h\t$n\t$r\t$a_het\t$a_theta\n" >> ../../data/fly/stats_popStat_table.txt
	printf "$t\t$c\t$G\t$h\t$n\t$r\t$mem\t$time\n" >> ../../data/fly/stats_compStat_table.txt


done
done
done
done
done
done
