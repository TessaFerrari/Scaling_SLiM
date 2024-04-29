# Script to make SLIM job script
# USAGE: ./make_slim_fly_main_coal.job.sh [n] [c] [G] [h] [rep]

cd /scratch1/tferrari/SlimBenchmark/fly

# Set n, the burn-in replicate number
n=${1}

# Set c, the scaling factor
c=${2}

# Set G, the genome size
G=${3}

# Set h, the deleterious dominance coefficient
h=${4}

# Set rep, the sim replicate number
rep=${5}

# Set t, the number of burn-in generations (grab length from burn-in stat file and multiply by c to unscale)
t=$((`awk 'BEGIN{FS="CYCLE "} /COALESCE/ {print $2+0}' ./out/burnCoal_scale${c}_gensize${G}/flyBench_burnin${n}_coal_c${c}_${G}.txt`*${c}))

# Make burn-in script
cat > ./scripts/temp/flyBench_coal_c${c}_${G}_h${h}_n${n}_rep${rep}.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 			// total chromosome length

	initializeMutationRate(6.2182e-9*${c});		// Deleterious mutation rate with a total mutation rate of 8.4e-9 and a ratio of neu to del of 1:2.85
	initializeMutationType("m2", ${h}, "g", -1.33e-4*${c}, 0.35); 	// Deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", m2, 1);			// Neutral will be overlaid later

	initializeGenomicElement(g1, 0, ${G}-1); 	// DFE from Huber et al. 2017 and recombination rate from Comeron et al. 2012
	initializeRecombinationRate(2.06e-8*${c});	// Demography from Sheehan and Song, 2016 - Three Epoch African Population
}

// Create the ancestral African population by reading burn-in trees file
1 early() { 
	
	defineConstant("simID", getSeed());
	sim.readFromPopulationFile("/scratch1/tferrari/SlimBenchmark/fly/out/burnCoal_scale${c}_gensize${G}/flyBench_burnin${n}_coal_c${c}_${G}.trees"); 
}

// Bottleneck 2.2M years ago
$((${t} / ${c} + 1)) early() { p1.setSubpopulationSize(asInteger(145300/${c})); }

// Population expansion 200k years ago
$(( ${t}/${c} + 2000000/${c} )) early() { p1.setSubpopulationSize(asInteger(544200/${c})); }

// Keep track of generation number in log file every 1000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} + 2200000/${c} )) late() {
        
	if (community.tick % 1000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/fly/logs/gen/flyBench_coal_c${c}_${G}_h${h}_n${n}_rep${rep}.gen", paste(sim.cycle));
        }
}

// After reaching present day, save tree sequence
$(( ${t}/${c} + 2200000/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/fly/out/burnCoal_scale${c}_gensize${G}/flyBench_coal_c${c}_${G}_h${h}_n${n}_rep${rep}.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: Coal");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
}

EOM
