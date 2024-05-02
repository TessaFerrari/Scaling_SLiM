# Script to make SLIM job script

# Set n, the burn-in replicate number
n=${1}

# Set t, the number of burn-in generations
t=${2}

# Set c, the scaling factor
c=${3}

# Set G, the genome size
G=${4}

# Set h, the deleterious dominance coefficient
h=${5}

# Set rep, the sim replicate number
rep=${6}

cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/fly

# Make burn-in script
cat > ./temp/flyBench_stats_${t}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.job << EOM

initialize() {
	
	defineConstant("L", ${G}); 			// total chromosome length

	initializeMutationRate(0); 			// No mutations, only collecting stats
	initializeMutationType("m1", 0.5, "f", 0.0);	// Neutral mutation
	initializeMutationType("m2", ${h}, "g", -1.33e-4*${c}, 0.35);   // Deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", c(m1,m2), c(1,2.85));        // Ratio of neu to del is 1:2.85


	initializeGenomicElement(g1, 0, ${G}-1); 
	initializeRecombinationRate(2.06e-8*${c});
}

// Read in tree and output stats
1 late() { 

	sim.readFromPopulationFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/fly/burn${t}_scale${c}_gensize${G}/flyBench_${t}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.trees");
	
	// Output VCFs
	p1_ids = p1.individuals.genomes;
        p1_ids.outputVCF(filePath="/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/fly/burn${t}_scale${c}_gensize${G}/flyBench_${t}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.african.vcf");	
	
	// Calculate Summary Stats
	p1_het = calcHeterozygosity(p1.genomes);	// African
	p1_theta = calcWattersonsTheta(p1.genomes); 

	// Output summary stats
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${t}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
	catn( "// ********** African heterozygousity: " + p1_het);
	catn( "// ********** African theta: " + p1_theta);

}

EOM
