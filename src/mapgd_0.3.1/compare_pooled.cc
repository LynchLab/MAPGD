#include "compare_pooled.h"

int compare_pooled(int argc, char *argv[])
{
	std::string infile="";
	std::string outfile="";

	env_t env;

	int s=0;
	map_double_t a=0.0;
	map_double_t EMLMIN=0.0001;

	std::vector <size_t> pop;

	env.setname("mapgd cp");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Micheal Lynch");

	env.setdescription("compares allele frequencies between pooled population genomic data.");

	env.optional_arg('m',"minerror", &EMLMIN, &arg_setfloat_t, "please provide an interger", "sets minimum");
	env.optional_arg('M',"MAX", &s, &arg_setint, "please provide an interger", "sets maximum");

	env.optional_arg('p',"populations", &pop, &arg_setvectorint, "please provide a list of comma seperated integers", "choose populations to compare. Populations should be specified by comma seperated\n\t\t\t\tintigers (e.g. 1,2,3) and cannot contain spaces (e.g. 1, 2, 3 is bad). A range can\n\t\t\t\tbe specified by using a hyphen (e.g. 1-10 for populations 1 through 10) and hypenated\n\t\t\t\tranges can be strung together (e.g. 1-10,15-16) so long as ranges do not overlap.");

	env.optional_arg('a',"alpha", &a, &arg_setfloat_t, "please provide a float", "only print sites where at least one population differes from the meta population mean with a p-value less than alpha.");

	env.optional_arg('i',"in", &infile, &arg_setstr, "please provide a valid inpuit file", "specifies input file (default datain.txt)");
	env.optional_arg('o',"out", &outfile, &arg_setstr, "please provide a valid name and location for output", "specifies output file (default stdout.txt) ");

//	env.flag('P',"allpopulations", &allpop, &flag_set, " ?? ", "compates all populations");
	env.flag('h',"help", &env, &flag_help, "an error occured while displaying the help mesaage", "prints this message");
	env.flag('v',"version", &env, &flag_version, "an error occured while displaying the version mesaage", "prints the program version");

	if ( parsargs(argc, argv, env) ) printUsage(env);
	if ( !env.required_set() ) printUsage(env);

	std::ostream *out;
	std::ofstream outFile;
	out=&std::cout;

	profile pro;
	allele_stat site;

	/* Open the input file. */

	if (infile.size()!=0) {
		pro.open(infile.c_str(), std::fstream::in);
		if (not (pro.is_open() ) ) printUsage(env);
	} else {
		pro.open(std::fstream::in);
	}

	/* Open the output file. */
	if (outfile.size()!=0) {
		outFile.open(outfile.c_str(), std::fstream::out);
		if (!outFile) printUsage(env);
		out=&outFile;
	}

	/* ************************************************************************************************************ */

	map_double_t pmaj;						/* fraction of putative major reads among the total of major and minor */
	map_double_t eml, pml, *pmlP;				/* maximum-likelihood estimates of the error rate and major-allele frequency */

	map_double_t *llhoodP, *llhoodS;				/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 
	map_double_t llhoodPS, llhoodSS;				/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 

	map_double_t *llstat;					/* difference in the log likelihoods */
	map_double_t maxll;						/* difference in the log likelihoods */

	lnmultinomial multi(4);	

	pro.maskall();

	if ( pop.size()==0 ) { 
		pop.clear();
		for (size_t x=0; x<pro.size(); ++x) pop.push_back(x);
	};

	for (size_t x=0; x<pop.size(); ++x) unmask(pro.get_locus().get_quartet(x));

	llhoodP=new map_double_t[pop.size()];
	llhoodS=new map_double_t[pop.size()];
	pmlP=new map_double_t[pop.size()];
	llstat=new map_double_t[pop.size()];

	/* 1) GET THE READS FOR ALL POPULATIONS, AND SET THE METAPOPULATION TO A AND B. */
	*out << "#id1\tid2\tmajor\tminor\t";
	for (size_t x=0; x<pop.size(); ++x) *out << "Freq\tdll\t";
	*out << "FreqMETA\tERROR" << std::endl;

	locus line;	
	while (pro.read(line)!=EOF ){
		line.maskall();		//maskall quartets
		line.unmask(pop);	//unmask quartets in pop
		line.mask_low_cov(1);	//mask low coverage

		/* 2) Identify the counts for the Putative major and minor alleles in the metapopulation. */ 

		//line=pro.getsite();
		line.sort();

		/* 3) CALCULATE THE LIKELIHOOD OF THE POOLED POPULATION DATA. */
	
		/* Calculate the ML estimates of the major pooled allele frequency and the error rate. */

		pmaj = (map_double_t) line.get_count(0) / (map_double_t) ( line.get_count(1) + line.get_count(0) );
		eml = 1.5 *(map_double_t) ( line.get_coverage() - line.get_count(0) - line.get_count(1) ) / ( (map_double_t) line.get_coverage() );
			if (eml < EMLMIN) { eml = EMLMIN;}
		pml = (pmaj * (1.0 - (2.0 * eml / 3.0) ) ) - (eml / 3.0);
		pml = pml / (1.0 - (4.0 * eml / 3.0) );

		if (pml > 1.0) { pml = 1.0; }
		if (pml < 0.0) { pml = 0.0; }
 
		/* Calculate the likelihoods under the reduced model assuming no variation between populations. */

		site.freq=pml;
		site.error=eml;
		site.major=line.get_index(0);
		site.minor=line.get_index(1);
		multi.set(&polymorphicmodel, site);
		for (size_t x=0; x<pop.size(); ++x) llhoodP[x]=multi.lnprob(line.get_quartet(pop[x]).base );

		/* 4) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF POPULATION SUBDIVISION. */ 

		for (size_t x=0; x<pop.size(); ++x) {
			pmaj = (map_double_t) line.get_count(pop[x],0) / (map_double_t) ( line.get_count(pop[x],1) + line.get_count(pop[x],0) );
			pmlP[x] = (pmaj * (1.0 - (2.0 * eml / 3.0) ) ) - (eml / 3.0);
			pmlP[x] = pmlP[x] / (1.0 - (4.0 * eml / 3.0) );
			if (pmlP[x] > 1.0) { pmlP[x] = 1.0; }
			if (pmlP[x] < 0.0) { pmlP[x] = 0.0; }

			site.freq=pmlP[x];
			multi.set(&polymorphicmodel, site);
                        llhoodP[x]=multi.lnprob(line.get_quartet(pop[x]).base );
			multi.set(&monomorphicmodel, site);
			llhoodS[x]=multi.lnprob(line.get_quartet(pop[x]).base ); 
		};
		
		/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */
	
		maxll=0;
		for (size_t x=0; x<pop.size(); ++x){
			if (!line.get_quartet(x).masked){
				llhoodSS=0;
				llhoodPS=0;
				for (size_t y=0; y<pop.size(); ++y){
					llhoodSS+=llhoodS[y];
					if (x!=y) llhoodPS+=llhoodS[y];
				else llhoodPS+=llhoodP[y];
				};
				llstat[x] = fabs(2.0 * (llhoodSS - llhoodPS) );
				if (llstat[x]>maxll) maxll=llstat[x];
			}
		};
		if (std::max(maxll, map_double_t(0) )>=a){
			*out << std::fixed << std::setprecision(7) << pro.getids(line) << '\t' << line.get_name(0) << '\t' << line.get_name(1) << '\t';
			for (size_t x=0; x<pop.size(); ++x){
				if ( line.get_coverage(pop[x])==0) *out << std::fixed << std::setprecision(7) << "NA" << '\t' << 0.0 << '\t';
				else *out << std::fixed << std::setprecision(7) << pmlP[x] <<'\t' << llstat[x] << '\t';
			};
			/* Significance at the 0.05, 0.01, 0.001 levels requires the statistic, with 1 degrees of freedom, to exceed 3.841, 6.635, and 10.827, respectively. */
			*out << pml << '\t' << eml << std::endl;
		};
	}
	delete [] llhoodP;
	delete [] llhoodS;
	delete [] llstat;
	delete [] pmlP;
	if (outFile.is_open()) outFile.close();
	exit(0);
};
