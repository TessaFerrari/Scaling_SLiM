# Script to make SLIM job script
# USAGE: make_slim_fly_main_sim.job.sh [n] [t] [c] [G] [h] [rep]

# Set n, the burn-in replicate number
n=${1}

# Set t, the number of burn-in generations (t = xNe; x=[5,10,20])
x=${2}
t=$(( ${x}*652700 ))

# Set c, the scaling factor
c=${3}

# Set G, the genome size
G=${4}

# Set h, the deleterious dominance coefficient
h=${5}

# Set rep, the sim replicate number
rep=${6}

cd /scratch1/tferrari/SlimBenchmark/fly

# Make burn-in script
cat > ./scripts/temp/flyBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 			// total chromosome length

	initializeMutationRate(8.4e-9*${c});		// Total mutation rate of 8.4e-9
	initializeMutationType("m1", 0.5, "f", 0.0);	// Neutral mutation type
	initializeMutationType("m2", ${h}, "g", -1.33e-4*${c}, 0.35); 	// Deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", c(m1,m2), c(1,2.85));	// Ratio of neu to del is 1:2.85

	initializeGenomicElement(g1, 0, ${G}-1); 	// DFE from Huber et al. 2017 and recombination rate from Comeron et al. 2012
	initializeRecombinationRate(2.06e-8*${c});	// Demography from Sheehan and Song, 2016 - Three Epoch African Population
}

// Create the ancestral African population by reading burn-in trees file
1 early() { 
	
	defineConstant("simID", getSeed());
	sim.readFromPopulationFile("/scratch1/tferrari/SlimBenchmark/fly/out/burn${x}_scale${c}_gensize${G}/flyBench_burnin${n}_${x}Ne_c${c}_${G}.trees"); }

// Bottleneck 2.2M years ago
$((${t} / ${c} + 1)) early() { p1.setSubpopulationSize(asInteger(145300/${c})); }

// Population expansion 200k years ago
$(( ${t}/${c} + 2000000/${c} )) early() { p1.setSubpopulationSize(asInteger(544200/${c})); }

// Keep track of generation number in log file every 1000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} + 2200000/${c} )) late() {
        
	if (community.tick % 1000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/fly/logs/gen/flyBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.gen", paste(sim.cycle));
        }
}

// After reaching present day, save tree sequence
$(( ${t}/${c} + 2200000/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/fly/out/burn${x}_scale${c}_gensize${G}/flyBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${x}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
}

EOM
