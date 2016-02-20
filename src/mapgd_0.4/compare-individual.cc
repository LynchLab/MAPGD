#include "compare-pooled.h"

int comparePooled(int argc, char *argv[])
{
	std::string infile="";
	std::string outfile="";

	Environment env;

	int s=0;
	float_t a=0.0;
	float_t EMLMIN=0.0001;

	std::vector <int> pop;

	env.set_name("mapgd cp");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman and Micheal Lynch");

	env.set_author("compares allele frequencies between pooled population genomic data.");

	env.optional_arg('m',"minerror", &EMLMIN, &arg_setfloat_t, "please provide an interger", "sets minimum");
	env.optional_arg('M',"MAX", &s, &arg_setint, "please provide an interger", "sets maximum");

	env.optional_arg('p',"populations", &pop, &arg_setvectorint, "please provide a list of comma seperated integers", "choose populations to compare. Populations should be specified by comma seperated\n\t\t\t\tintigers (e.g. 1,2,3) and cannot contain spaces (e.g. 1, 2, 3 is bad). A range can\n\t\t\t\tbe specified by using a hyphen (e.g. 1-10 for populations 1 through 10) and hypenated\n\t\t\t\tranges can be strung together (e.g. 1-10,15-16) so long as ranges do not overlap.");

	env.optional_arg('a',"alpha", &a, &arg_setfloat_t, "please provide a float", "only print sites where at least one population differes from the meta population mean with a p-value less than alpha.");

	env.optional_arg('i',"in", &infile, &arg_setstr, "please provide a valid inpuit file", "specifies input file (default datain.txt)");
	env.optional_arg('o',"out", &outfile, &arg_setstr, "please provide a valid name and location for output", "specifies output file (default stdout.txt) ");

//	env.flag('P',"allpopulations", &allpop, &flag_set, " ?? ", "compates all populations");
	env.flag('h',"help", &env, &flag_help, "an error occured while displaying the help mesaage", "prints this message");
	env.flag('v',"version", &env, &flag_version, "an error occured while displaying the version mesaage", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env);
	if ( !env.required_set() ) print_usage(env);

	std::ostream *out;
	std::ofstream outFile;
	out=&std::cout;

	profile pro;
	Allele site;

	/* Open the input file. */

	if (infile.size()!=0) {if (pro.open(infile.c_str(), "r")==NULL) {print_usage(env);} }
	else pro.open("r");

	/* Open the output file. */
	if (outfile.size()!=0) {
		outFile.open(outfile, std::ofstream::out);
		if (!outFile) print_usage(env);
		out=&outFile;
	}

	/* ************************************************************************************************************ */

	float_t pmaj;						/* fraction of putative major reads among the total of major and minor */
	float_t eml, pml, *pmlP;				/* maximum-likelihood estimates of the error rate and major-allele frequency */

	float_t *llhoodP, *llhoodS;				/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 
	float_t llhoodPS, llhoodSS;				/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 

	float_t *llstat;					/* difference in the log likelihoods */
	float_t maxll;						/* difference in the log likelihoods */

	lnmultinomial multi(4);	

	pro.maskall();

	if ( pop.size()==0 ) { 
		pop.clear();
		for (size_t x=0; x<pro.size(); ++x) pop.push_back(x);
	};

	for (size_t x=0; x<pop.size(); ++x) pro.unmask(pop[x]);

	llhoodP=new float_t[pop.size()];
	llhoodS=new float_t[pop.size()];
	pmlP=new float_t[pop.size()];
	llstat=new float_t[pop.size()];

	/* 1) GET THE READS FOR ALL POPULATIONS, AND SET THE METAPOPULATION TO A AND B. */
	*out << "#id1\tid2\tmajor\tminor\t";
	for (size_t x=0; x<pop.size(); ++x) *out << "Freq\tdll\t";
	*out << "FreqMETA\tERROR" << std::endl;

	Locus line;	
	while (pro.read(line)!=EOF ){
		for (size_t x=0; x<pop.size(); ++x) pro.unmask(pop[x]);
		for (size_t x=0; x<pop.size(); ++x) line.unmask(pop[x]);
		for (size_t x=0; x<pop.size(); ++x) if (line.getcoverage(x)==0) line.sample[x].masked=true;

		/* 2) Identify the counts for the Putative major and minor alleles in the metapopulation. */ 

		//line=pro.getsite();
		line.sort();

		/* 3) CALCULATE THE LIKELIHOOD OF THE POOLED POPULATION DATA. */
	
		/* Calculate the ML estimates of the major pooled allele frequency and the error rate. */

		pmaj = (float_t) line.getcount(0) / (float_t) ( line.getcount(1) + line.getcount(0) );
		eml = 1.5 *(float_t) ( line.getcoverage() - line.getcount(0) - line.getcount(1) ) / ( (float_t) line.getcoverage() );
			if (eml < EMLMIN) { eml = EMLMIN;}
		pml = (pmaj * (1.0 - (2.0 * eml / 3.0) ) ) - (eml / 3.0);
		pml = pml / (1.0 - (4.0 * eml / 3.0) );

		if (pml > 1.0) { pml = 1.0; }
		if (pml < 0.0) { pml = 0.0; }
 
		/* Calculate the likelihoods under the reduced model assuming no variation between populations. */

		site.freq=pml;
		site.error=eml;
		site.major=line.getindex(0);
		site.minor=line.getindex(1);
		multi.set(&polymorphicmodel, site);
		for (size_t x=0; x<pop.size(); ++x) llhoodP[x]=multi.lnprob(line.get_quartet(pop[x]).base );

		/* 4) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF POPULATION SUBDIVISION. */ 

		for (size_t x=0; x<pop.size(); ++x) {
			pmaj = (float_t) line.getcount(pop[x],0) / (float_t) ( line.getcount(pop[x],1) + line.getcount(pop[x],0) );
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
			if (!line.sample[x].masked){
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
		if (std::max(maxll, float_t(0) )>=a){
			*out << std::fixed << std::setprecision(7) << pro.getids(line) << '\t' << line.getname(0) << '\t' << line.getname(1) << '\t';
			for (size_t x=0; x<pop.size(); ++x){
				if ( line.getcoverage(pop[x])==0) *out << std::fixed << std::setprecision(7) << "NA" << '\t' << 0.0 << '\t';
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
