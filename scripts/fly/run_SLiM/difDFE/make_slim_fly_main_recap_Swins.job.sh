# Script to make SLIM job script
# USAGE: ./make_slim_fly_main_recap.job.sh [n] [c] [G] [h] [rep]

# Set n, the burn-in replicate number
n=${1}

# Set t to 0 since sim will be recapitated
t=0

# Set c, the scaling factor
c=${2}

# Set G, the genome size
G=${3}

# Set h, the deleterious dominance coefficient
h=${4}

# Set rep, the sim replicate number
rep=${5}

cd /scratch1/tferrari/SlimBenchmark/fly

# Make burn-in script
cat > ./scripts/temp/flyBench_recap_c${c}_${G}_h${h}_n${n}_rep${rep}_Swins.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 			// total chromosome length

	// Create evenly distributed 1kb windows of selection in the genome, so that 1/5 of the genome undergoes selection
	ends = sort(c(seq(3999, L-1, by=5000),seq(4999, L-1, by=5000)));
	rates = rep(c(0,6.2182e-9*${c}), asInteger(L/5000) );
	initializeMutationRate(rates, ends);		// Deleterious mutation rate with a total mutation rate of 8.4e-9 and a ratio of neu to del of 1:2.85
	
	initializeMutationType("m2", ${h}, "g", -1.33e-4*${c}, 0.35); 	// Deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", m2, 1);			// Neutral will be overlaid later

	initializeGenomicElement(g1, 0, ${G}-1); 	// DFE from Huber et al. 2017 and recombination rate from Comeron et al. 2012
	initializeRecombinationRate(2.06e-8*${c});	// Demography from Sheehan and Song, 2016 - Three Epoch African Population

}

// Create the ancestral African population by reading burn-in trees file
1 early() { 

	defineConstant("simID", getSeed());
	sim.addSubpop("p1", asInteger(652700/${c})); 
}

// Bottleneck 2.2M years ago
$((${t} / ${c} + 1)) early() { p1.setSubpopulationSize(asInteger(145300/${c})); }

// Population expansion 200k years ago
$(( ${t}/${c} + 2000000/${c} )) early() { p1.setSubpopulationSize(asInteger(544200/${c})); }

// Keep track of generation number in log file every 10000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} + 2200000/${c} )) late() {
        
	if (community.tick % 10000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/fly/logs/gen/flyBench_recap_c${c}_${G}_h${h}_n${n}_rep${rep}_Swins.gen", paste(sim.cycle));
        }
}

// After reaching present day, save tree sequence
$(( ${t}/${c} + 2200000/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/fly/out/burnRecap_scale${c}_gensize${G}/flyBench_recap_c${c}_${G}_h${h}_n${n}_rep${rep}_Swins.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: Recap");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
}

EOM
