# Script to make SLIM job script

# Set n, the burn-in replicate number
n=${1}

# Set t, the number of burn-in generations (t = xNe; x=[5,10,20])
x=${2}
t=$(( ${x}*7310 )) # unscaled

# Set c, the scaling factor
c=${3}

# Set G, the genome size
G=${4}

cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/human

# Make burn-in script
cat > ./temp/humanBench_burnin${n}_${x}Ne_c${c}_${G}.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 			// Total chromosome length

	initializeMutationRate(2.124e-8*${c});		// Nine tenths of the original mutation rate 2.36e-8 (neutral only)
	initializeMutationType("m1", 0.5, "f", 0.0);	// Neutral mutations
	initializeGenomicElementType("g1", m1, 1);

	initializeGenomicElement(g1, 0, ${G}-1); 
	initializeRecombinationRate(1e-8*${c});		// Recombination rate from SLiM Recipe 5.4 - The Gravel et al. (2011) model of human evolution
}

// Create the ancestral African population
1 early() { 
	
	defineConstant("simID", getSeed());
	sim.addSubpop("p1", asInteger(7310/${c}) ); 
}

// Keep track of generation number in log file every 1000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} )) late() {
        
	if (community.tick % 1000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/gen_logs/human/humanBench_burnin${n}_${x}Ne_c${c}_${G}.gen", paste(sim.cycle));
        }
}

// After burn-in, save tree sequence
$(( ${t}/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${x}_scale${c}_gensize${G}/humanBench_burnin${n}_${x}Ne_c${c}_${G}.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${x}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
}

EOM
