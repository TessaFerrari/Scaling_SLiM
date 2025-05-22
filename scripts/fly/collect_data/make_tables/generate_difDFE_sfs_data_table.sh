# Set number of replicates per burn-in
rep=1 

cd /scratch1/tferrari/SlimBenchmark/fly/out

out=../data/difDFE_subsamp50_sfs_table.txt

#Print sfs header
printf "Population\tBurnInType\tScalingFactor\tGenomeSize\tDomCoefficent\tBurnNum\tSimNum\tDFEtype\tTotNumMuts" > ${out}
for i in {1..50}; do printf "\t$i"; done >> ${out}
printf "\n" >> ${out}

for pop in african; do
for t in Recap; do
if [ "$t" != "Coal" ] && [ "$t" != "Recap" ]; then x=$t"Ne"; elif [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

for c in 1000 500 100; do
for G in 1e5 1e6 1e7; do
#if [[ "$G" == "1e5" ]]; then z=$(($rep*10)); else z=$rep; fi

for h in 0.5; do
for n in 1; do
for DFEtype in lowS Swins; do

	file=burn${t}_scale${c}_gensize${G}/sfs/flyBench_${x}_c${c}_${G}_h${h}_n${n}_${DFEtype}.${pop}.txt

	# Print SFS without header
	cat $file | tail -n +2 >> ${out}

done
done
done
done
done
done
done

