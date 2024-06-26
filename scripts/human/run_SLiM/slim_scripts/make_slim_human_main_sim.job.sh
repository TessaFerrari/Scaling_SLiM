# Script to make SLIM job script

# Set n, the burn-in replicate number
n=${1}

# Set t, the number of burn-in generations (t = x*Ne; x=[5,10,20])
x=${2}
t=$(( ${x}*7310 ))

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
cat > ./temp/humanBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.job << EOM

initialize() {
	
	initializeTreeSeq();
	defineConstant("L", ${G}); 					// Total chromosome length

	initializeMutationRate(2.36e-8*${c});				// Full mutation rate (for both neu and del mutations)
	initializeMutationType("m1", 0.5, "f", 0.0);			// Neutral mutation
	initializeMutationType("m2", ${h}, "g", -0.03*${c}, 0.2); 	// Mild deleterious mutation (recessive or additive)
	initializeGenomicElementType("g1", c(m1,m2), c(9,1));		// Ratio of Neu:Del is 9:1

	initializeGenomicElement(g1, 0, ${G}-1); 
	initializeRecombinationRate(1e-8*${c});				// Recombination rate from SLiM Recipe 5.4 - The Gravel et al. (2011) model of human evolution
}

// Create the ancestral African population by reading burn-in trees file
1 late() { 
	defineConstant("simID", getSeed());
	sim.readFromPopulationFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${x}_scale${c}_gensize${G}/humanBench_burnin${n}_${x}Ne_c${c}_${G}.trees"); }

// Expand the African population to 14474
// This occurs 148000 years (5920) generations ago
$((${t} / ${c} + 1)) early() { p1.setSubpopulationSize(asInteger(14474/${c})); }

// Split non-Africans from Africans and set up migration between them
// This occurs 51000 years (2040 generations) ago
$(( ${t}/${c} + 3880/${c} )) early() {
	sim.addSubpopSplit("p2", asInteger(1861/${c}), p1);
	p1.setMigrationRates(c(p2), c(15e-5));
	p2.setMigrationRates(c(p1), c(15e-5));
}

// Split p2 into European and East Asian subpopulations
// This occurs 23000 years (920 generations) ago
$(( ${t}/${c} + 5000/${c} )) early() {
	sim.addSubpopSplit("p3", asInteger(554/${c}), p2);
	p2.setSubpopulationSize(asInteger(1032/${c}));  // reduce European size

	// Set migration rates for the rest of the simulation
	p1.setMigrationRates(c(p2, p3), c(2.5e-5, 0.78e-5));
	p2.setMigrationRates(c(p1, p3), c(2.5e-5, 3.11e-5));
	p3.setMigrationRates(c(p1, p2), c(0.78e-5, 3.11e-5));
}

// Set up exponential growth in Europe and East Asia
// Where N(0) is the base subpopulation size and t = gen - 57080:
//    N(Europe) should be int(round(N(0) * e^(0.0038*t)))
//    N(East Asia) should be int(round(N(0) * e^(0.0048*t)))
$(( ${t}/${c} + 5000/${c} )):$(( ${t}/${c} + 5920/${c} )) early() {
	t = community.tick - $(( ${t}/${c} + 5000/${c} ));
	p2_size = round(1032 * exp(0.0038 * t * ${c}));
	p3_size = round(554 * exp(0.0048 * t * ${c}));
	
	p2.setSubpopulationSize(asInteger(p2_size/${c}));
	p3.setSubpopulationSize(asInteger(p3_size/${c}));
}

// Keep track of generation number in log file every 1000 generations until burn-in end (overwrites itself)
1:$(( ${t}/${c} + 5920/${c} )) late() {
        
	if (community.tick % 1000 == 0){
                writeFile("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/gen_logs/human/humanBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.gen", paste(sim.cycle));
        }
}

// After reaching present day, save tree sequence
$(( ${t}/${c} + 5920/${c} )) late() {
	
	sim.treeSeqOutput("/scratch1/tferrari/SlimBenchmark/Scaling_SLiM/out/human/burn${x}_scale${c}_gensize${G}/humanBench_${x}Ne_c${c}_${G}_h${h}_n${n}_rep${rep}.trees");
	catn( "// ********** Initial random number seed: " + simID);
	catn( "// ********** Burn-in replicate number: ${n}");
	catn( "// ********** Burn-in type: ${x}");
	catn( "// ********** Scaling factor: ${c}");
	catn( "// ********** Genome size: ${G}");
	catn( "// ********** Dominance coefficient: ${h}");
	catn( "// ********** Simulation replicate number: ${rep}");
}

EOM
