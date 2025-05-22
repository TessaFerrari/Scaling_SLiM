# Set bins, the number of bins
nbins=10

# Set dist, the maximum distance between snps
dist=1000

cd /scratch1/tferrari/SlimBenchmark/fly/out

out=../data/dist${dist}_${nbins}bins_ld_table.txt

#Print ld header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum" > ${out}
for (( i=1; i <= $nbins; i++ )); do printf "\tb"$i"_nSNPs\tb"$i"_sumR2\tb"$i"_sumD\tb"$i"_Dp"; done >> ${out}
printf "\n" >> ${out}


for pop in african; do
for t in 5 10 20 Coal Recap; do
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

for c in 1000 500 100; do
for G in 1e5 1e6 1e7; do

for h in 0.5; do
for n in 1 2; do

	file=burn${t}_scale${c}_gensize${G}/ld/dist${dist}_${nbins}bins_flyBench_${x}_c${c}_${G}_h${h}_n${n}.${pop}.txt

	# Print SFS without header
	cat $file | tail -n +2 >> ${out}

done
done
done
done
done
done
