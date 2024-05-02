# Script to make SLIM job script

# Set n, the burn-in replicate number
n=${1}

# Set t, the tree type (Coal or Recap)
t=${2}
if [ "$t" == "Coal" ]; then x=coal; else x=recap; fi

# Set c, the scaling factor
c=${3}

# Set G, the genome size
G=${4}

# Set h, the deleterious dominance coefficient
h=${5}

# Set rep, the sim replicate number
rep=${6}

cd /scratch1/tferrari/SlimBenchmark/Scaling_SLiM/scripts/human

# Make burn-in script
cat > ./temp/humanBench_stats_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.job << EOM

initialize() {
	
	defineConstant("L", ${G}); 			// total chromosome length

	initializeMutationRate(0); 			// No mutations, only collecting stats
	initializeMutationType("m1", 0.5, "f", 0.0);	// Neutral mutation
	initializeMutationType("m2", ${h}, "g", -0.03*${c}, 0.2); 	// Mild deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", c(m1,m2), c(9,1));

	initializeGenomicElement(g1, 0, ${G}-1); 
	initializeRecombinationRate(1e-8*${c});		// recombination rate from SLiM Recipe 5.4 - The Gravel et al. (2011) model of human evolution
}

// Read in tree and output stats
1 late() { 

	sim.readFromPopulationFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${t}_scale${c}_gensize${G}/humanBench_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.over.trees");
	
	// Output VCFs
	p1_ids = p1.individuals.genomes;
        p1_ids.outputVCF(filePath="/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${t}_scale${c}_gensize${G}/humanBench_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.african.vcf");	
	p2_ids = p2.individuals.genomes;
        p2_ids.outputVCF(filePath="/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${t}_scale${c}_gensize${G}/humanBench_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.european.vcf");
	p3_ids = p3.individuals.genomes;
        p3_ids.outputVCF(filePath="/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${t}_scale${c}_gensize${G}/humanBench_${x}_c${c}_${G}_h${h}_n${n}_rep${rep}.eastasian.vcf");

	// Calculate Summary Stats
	p1_het = calcHeterozygosity(p1.genomes);	// African
	p1_theta = calcWattersonsTheta(p1.genomes); 
        p2_het = calcHeterozygosity(p2.genomes);	// European
        p2_theta = calcWattersonsTheta(p2.genomes);
        p3_het = calcHeterozygosity(p3.genomes);	// East Asian
        p3_theta = calcWattersonsTheta(p3.genomes);

	// Output summary stats
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${t}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
	catn( "// ********** African heterozygousity: " + p1_het);
	catn( "// ********** African theta: " + p1_theta);
	catn( "// ********** European heterozygousity: " + p2_het);
	catn( "// ********** European theta: " + p2_theta);
	catn( "// ********** East Asian heterozygousity: " + p3_het);
	catn( "// ********** East Asian theta: " + p3_theta);

}

EOM
