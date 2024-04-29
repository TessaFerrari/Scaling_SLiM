# Script to make SLIM job script
# USAGE: make_slim_fly_burnin.job.sh [n] [t] [c] [G]

# Set n, the burn-in replicate number
n=${1}

# Set t, the number of burn-in generations (t = xNe; x=[5,10,20])
x=${2}
t=$(( ${x}*652700 ))

# Set c, the scaling factor
c=${3}

# Set G, the genome size
G=${4}

cd /scratch1/tferrari/SlimBenchmark/fly

# Make burn-in script
cat > ./scripts/temp/flyBench_burnin${n}_${x}Ne_c${c}_${G}.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 			// total chromosome length

	initializeMutationRate(2.1818e-9*${c});		// Neutral mutation rate with an original mutation rate of 8.4e-9 and a neu to del ratio of 1:2.85
	initializeMutationType("m1", 0.5, "f", 0.0);	// Neutral mutation type
	initializeGenomicElementType("g1", m1, 1);

	initializeGenomicElement(g1, 0, ${G}-1); 	// DFE from Huber et al. 2017 and recombination rate from Comeron et al. 2012
	initializeRecombinationRate(2.06e-8*${c});	// Demography from Sheehan and Song, 2016 - Three Epoch African Population
}

// Create the ancestral African population
1 early() { 
	
	defineConstant("simID", getSeed());
	sim.addSubpop("p1", asInteger(652700/${c}) ); 
}

// Keep track of generation number in log file every 10000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} )) late() {
        
	if (community.tick % 10000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/fly/logs/gen/flyBench_burnin${n}_${x}Ne_c${c}_${G}.gen", paste(sim.cycle));
        }
}

// After burn-in, save tree sequence
$(( ${t}/${c} )) early() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/fly/out/burn${x}_scale${c}_gensize${G}/flyBench_burnin${n}_${x}Ne_c${c}_${G}.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${x}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
}

EOM
