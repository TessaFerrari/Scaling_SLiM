# Set number of replicates per burn-in
rep=5 

cd /scratch1/tferrari/SlimBenchmark/fly/out

#out=../data/mutStats_table.txt
out=../data/mutStatsStart_table.txt

#Print header
#printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tTot.MeanS\tDel.MeanS\tTot.MeanTO\tDel.MeanTO\tNeu.MeanTO\tCurrentTick\n" > ${out}
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tStartingTick\n" > ${out}

for pop in african; do
for t in 5 10 20 Coal Recap; do
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

for c in 1000 500 100; do
for G in 1e5 1e6 1e7; do
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

for h in 0.5; do
for n in 1 2; do

	#file=burn${t}_scale${c}_gensize${G}/mutStats/flyBench_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt
	file=burn${t}_scale${c}_gensize${G}/mutStats/flyBench_start_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

	# Print without header
	cat $file | tail -n +2 >> ${out}

done
done
done
done
done
done
