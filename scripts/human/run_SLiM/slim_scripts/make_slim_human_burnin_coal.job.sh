# Script to make SLIM job script

# Set n, the burn-in replicate number
n=${1}

# Set c, the scaling factor
c=${2}

# Set G, the genome size
G=${3}

# Set t, the unscaled maximum burn-in length (30Ne)
t=$((30*7310))

cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/human

# Make burn-in script
cat > ./temp/humanBench_burnin${n}_coal_c${c}_${G}.job << EOM

initialize() {
	
	initializeTreeSeq(checkCoalescence=T);
	defineConstant("L", ${G}); 			// Total chromosome length

	initializeMutationRate(0);			// Set mutation rate to 0, neutral mutations will be overlaid later
	initializeMutationType("m1", 0.5, "f", 0.0);
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
                writeFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/gen_logs/human/humanBench_burnin${n}_coal_c${c}_${G}.gen", paste(sim.cycle));
        }
}

// Once tree has coalesced, save tree sequence and end simulation
1:$(( ${t}/${c} )) late() {
	if (sim.treeSeqCoalesced()) {
		sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burnCoal_scale${c}_gensize${G}/humanBench_burnin${n}_coal_c${c}_${G}.trees");
		catn("COALESCED AT CYCLE " + sim.cycle + " (" + sim.cycle*${c}/7310 + "Ne)");
		catn( "// ********** Initial random number seed: " + simID);
		catn( "// ********** Burn-in replicate number: ${n}");
		catn( "// ********** Burn-in type: Coal");
		catn( "// ********** Scaling factor: ${c}");
		catn( "// ********** Genome size: ${G}");
		sim.simulationFinished();
	}
}

// If tree hasn't coalesced after 30Ne generations, save tree sequence
$(( ${t}/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burnCoal_scale${c}_gensize${G}/humanBench_burnin${n}_coal_c${c}_${G}.trees");
	catn("NO COALESCENCE BY CYCLE $(( ${t}/${c} )) (30Ne)");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: Coal");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
}

EOM
