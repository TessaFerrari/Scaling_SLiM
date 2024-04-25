# Set number of replicates per burn-in
rep=5 

cd /scratch1/tferrari/SlimBenchmark/fly/out

printf "BurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tADel\tANeu\tEDel\tENeu\tEADel\tEANeu\n" > ../data/mutCounts_table.txt

for t in 5 10 20 Coal Recap; do
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

for c in 1000 500 100; do
for G in 1e5 1e6 1e7; do
if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

for h in 0.5; do
for n in 1 2; do

for (( r=1; r <= $z; r++ )); do

	file=burn${t}_scale${c}_gensize${G}/mut_counts/flyBench_${x}_c${c}_${G}_h${h}_n${n}_rep${r}.mutcounts.txt

	# Mutation counts
	a_del=`awk -F": " '/african deleterious:/ {print $2}' $file`
	a_neu=`awk -F": " '/african neutral:/ {print $2}' $file`

	printf "$t\t$c\t$G\t$h\t$n\t$r\t$a_del\t$a_neu\n" >> ../data/mutCounts_table.txt

done
done
done
done
done
done
