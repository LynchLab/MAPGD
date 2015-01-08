#include 	"interface.h"
#include 	"comparePooled.h"
#include 	<iostream>

#define EMLMIN	0.00001
#define LNFMAX	10000

float_t lnf[LNFMAX];

void lnfactinit(void){
	lnf[0]=0;
	for (int x=1; x<LNFMAX; ++x) lnf[x]=lnf[x-1]+log(x);
}

float_t lnfact(count_t x){

	if (x<LNFMAX) return lnf[x];
	else{
		float_t fx=x;
		return  fx * log( fx - 1.0) + log(pow(2.0 * M_PI * fx, 0.5 )  ); 
	} 

}

int comparePooled(int argc, char *argv[])
{
	std::string infile="datain.txt";
	std::string outfile="dataout.txt";

	env_t env;
	bool sane=false;

	int t=0;
	int s=0;
	int pop[2];

	env.setname("mapgd cp");
	env.setver("1.0");
	env.setauthor("Micheal Lynch");

	env.setdescription("compares allele frequencies between pooled population genomic data");

	env.optional_arg('m',"min", &t, &arg_setint, "please provide an interger", "sets minimum");
	env.optional_arg('M',"MAX", &s, &arg_setint, "please provide an interger", "sets maximum");

	env.optional_arg('p',"population", &pop, &arg_set2int, "please provide two intergers", "choose two populations to compare");

//	env.optional_arg('m',"mgd", &sitefile, &arg_setstr, "please provide a valid mgd file", "a file that provides likelihood of polymorphism in each population.");

	env.optional_arg('i',"in", &infile, &arg_setstr, "please provide a valid inpuit file", "specifies input file (default datain.txt)");
	env.optional_arg('o',"out", &outfile, &arg_setstr, "please provide a valid name and location for output", "specifies output file (default stdout.txt) ");

	env.flag('h',"help", &env, &flag_help, "an error occured while displaying the help mesaage", "prints this message");
	env.flag('s',"sane", &sane, &flag_set, "takes no argument", "set default in/out to the stdin and stdout.");
	env.flag('v',"version", &env, &flag_version, "an error occured while displaying the version mesaage", "prints the program version");

	if ( parsargs(argc, argv, env) ) printUsage(env);
	if ( !env.required_set() ) printUsage(env);

	if ( sane ) { infile=""; outfile=""; }
	
	lnfactinit();

	std::ostream *out;
	out=&std::cout;
	profile pro;

	/* Open the input file. */

	std::cout << infile.size() << std::endl;
	if (infile.size()!=0) {if (pro.open(infile.c_str(), 'r')==NULL) {printUsage(env);} }
	else pro.open('r');
/*
	if (outfile.size()!=0) {
		outFile.open(outfile, 'w');
		if (!outFile) printUsage(env);
		out=&outFile;
	}*/


	/* Open the input file. */


	/* ************************************************************************************************************ */

	double pmaj;						/* fraction of putative major reads among the total of major and minor */
	double pmajA, pmajB;

	double eml, pml;					/* maximum-likelihood estimates of the error rate and major-allele frequency */
	double pmlA, pmlB;

	double totests;						/* sums to use in estimating mean and standard deviation of the estimates over all runs */

	double emltheta;

	double lnterm;
	float_t lnpml, ln1pml;

	float_t exphimaj[20001];					/* terms used in the likelihood test statistic */
	float_t exphimin[20001];

	float_t lhoodP1, lhoodS1;				/* likelihoods of the data in population 1 under the assumptions of pooled or distinct allele frequencies */
	float_t lhoodP2, lhoodS2;				/* likelihoods of the data in population 2 under the assumptions of pooled or distinct allele frequencies */

	float_t llhoodP, llhoodS;				/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 

	count_t minig, maxig, ig;					/* limits to allele frequency count used in the likelihood analysis, to streamline the computational time */

	count_t coverage;

	double llstat;						/* difference in the log likelihoods */

	count_t a=pop[0]-1, b=pop[1]-1;

	/* Precalculate the binomial terms for the full likelihood. */

	
	/* ********** Set the minor-allele frequencies for the first population at which the computations will be done. ************* */
	/* *alculate the starting and ending points for the allele counts used in the likelihood function. */

	while (pro.read()!=EOF ){

		/* 1) GET THE READS FOR ALL POPULATIONS, AND SET THE METAPOPULATION TO A AND B. */

		pro.maskall();
		pro.unmask(a);
		pro.unmask(b);

		coverage=pro.getcoverage();
		minig = 0;
		maxig=coverage;

		/* 4) IDENTIFY THE OVERALL MAJOR AND MINOR ALLELES. */ 

		/* Identify the counts for the Putative major and minor alleles	for the pooled pair of populations. */ 

		pro.sort();

		/* 5) CALCULATE THE LIKELIHOOD OF THE JOINTLY POOLED POPULATIONS. */
	
		/* Calculate the ML estimates of the major pooled allele frequency and the error rate. */

		pmaj = (float_t) pro.getcount(0) / (float_t) ( pro.getcount(1) + pro.getcount(0) );
		eml = 1.5 *(float_t) ( coverage - pro.getcount(0) - pro.getcount(1) ) / ( (float_t) coverage );
			if (eml == 0.0) { eml = EMLMIN;}
		pml = (pmaj * (1.0 - (2.0 * eml / 3.0) ) ) - (eml / 3.0);
		pml = pml / (1.0 - (4.0 * eml / 3.0) );
		emltheta = 1.0 - (4.0 * eml / 3.0);

		if (pml >= 1.0) { pml = 0.999999999; }

//		std::cout << EMLMIN << ", " << coverage << ", " << pro.getcount(0) << ", " << pro.getcount(1) << std::endl;	
//		std::cout << pmaj << ", :: " << eml << ", " << pml << ", " << emltheta << std::endl;
	
		/* Calculate the likelihoods under the reduced model assuming no variation between populations. */
			
		/* Calculate the component for the first population. */

		count_t coverageA=pro.getcoverage(a);
		count_t coverageB=pro.getcoverage(b);

		lhoodP1 = 0.0;
		
		for (ig = minig; ig <= coverageA; ++ig) {
			exphimaj[ig] = log(((((double) ig) / coverageA) * emltheta) + (eml / 3.0));
			exphimin[ig] = log(1.0 - eml - ((((double) ig) / coverageA) * emltheta)); }

		if ( (pml *(1.0-pml)) > 0.0) {
		
			lnpml = log(pml);
			ln1pml = log(1.0 - pml);

			for (ig = minig; ig <= coverageA; ++ig) {
				lnterm = lnfact(coverageA) - lnfact(ig) - lnfact(coverageA-ig) + (((double) ig) * lnpml) + (((double) (coverageA-ig)) * ln1pml)
					+ (((double) pro.getcount(a,0)) * exphimaj[ig]) +  (((double) pro.getcount(a,1)) * exphimin[ig]); 
				lhoodP1 = lhoodP1 + exp(lnterm); }
		}
		
		for (ig = minig; ig <= coverageB; ++ig) {
			exphimaj[ig] = log(((((double) ig) / coverageB) * emltheta) + (eml / 3.0));
			exphimin[ig] = log(1.0 - eml - ((((double) ig) / coverageB) * emltheta)); }

		/* Calculate the component for the second population. */

		lhoodP2 = 0.0;
		if ( (pml *(1.0-pml)) > 0.0) {
			for (ig = minig; ig <= coverageB; ++ig) {
				lnterm = lnfact(coverageB) - lnfact(ig) - lnfact(coverageB-ig) + (((double) ig) * lnpml) + (((double) (coverageB-ig)) * ln1pml)
					+ (((double) pro.getcount(b,0)) * exphimaj[ig]) +  (((double) pro.getcount(b,1)) * exphimin[ig]); 
				lhoodP2 = lhoodP2 + exp(lnterm); }
		}
		
		if ( (pml *(1.0-pml)) <= 0.0) {
			llhoodP = 0.0;}
		else { llhoodP = log(lhoodP1) + log(lhoodP2); }

		/* 6) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF POPULATION SUBDIVISION. */ 

		pmajA = ((double) pro.getcount(a,0)) / ( ((double) pro.getcount(a,0)) + ((double) pro.getcount(a,1)) );
		pmajB = ((double) pro.getcount(b,0)) / ( ((double) pro.getcount(b,0)) + ((double) pro.getcount(b,1)) );
	
		pmlA = (pmajA * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
		pmlA = pmlA / (1.0 - (4.0 * eml / 3.0));
		pmlB = (pmajB * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
		pmlB = pmlB / (1.0 - (4.0 * eml / 3.0)); 

		if (pmlA >= 1.0) { pmlA = 0.999999999; }
		if (pmlB >= 1.0) { pmlB = 0.999999999; }
		
		//*out << pmlA << ", " << pmlB << std::endl;

		/* Calculate the component for the first population. */

		lhoodS1 = 0.0;
		lnpml = log(pmlA);
		ln1pml = log(1.0 - pmlA);

		for (ig = minig; ig <= coverageA; ++ig) {
			exphimaj[ig] = log(((((double) ig) / coverageA) * emltheta) + (eml / 3.0));
			exphimin[ig] = log(1.0 - eml - ((((double) ig) / coverageA) * emltheta)); }
	
		if (  ( (pmlA * (1.0 - pmlA) ) != 0.0)  ) {
			for (ig = minig; ig <= coverageA; ++ig) {
				//*out << "lhood:" << lhoodS1 << std::endl;
				lnterm = lnfact(coverageA) - lnfact(ig) - lnfact(coverageA-ig) + ( ( (double) ig) * lnpml) + ( ( (double) (coverageA-ig)) * ln1pml)
					+ ( ((double) pro.getcount(a,0)) * exphimaj[ig]) +  (((double) pro.getcount(a,1)) * exphimin[ig]); 
				lhoodS1 = lhoodS1 + exp(lnterm); } } 
				//*out << "done:" << lhoodS1 << std::endl;
		
		/* Calculate the component for the second population. */

		lhoodS2 = 0.0;
		lnpml = log(pmlB);
		ln1pml = log(1.0 - pmlB);

		for (ig = minig; ig <= coverageB; ++ig) {
			exphimaj[ig] = log(((((double) ig) / coverageB) * emltheta) + (eml / 3.0));
			exphimin[ig] = log(1.0 - eml - ((((double) ig) / coverageB) * emltheta)); }
		
		if ( ( (pmlB * (1.0 - pmlB) ) != 0.0) ) { 
			for (ig = minig; ig <= coverageB; ++ig) {
				//*out << "lhood:" << lhoodS2 << std::endl;
				lnterm = lnfact(coverageB) - lnfact(ig) - lnfact(coverageB-ig) + (((double) ig) * lnpml) + (((double) (coverageB-ig)) * ln1pml)
					+ ( ((double) pro.getcount(b,0)) * exphimaj[ig]) +  (((double) pro.getcount(b,1)) * exphimin[ig]); 
				lhoodS2 = lhoodS2 + exp(lnterm); } }

		if ( ( (lhoodS1 * lhoodS2) == 0.0) ) {
			llhoodS = 0.0;
		} else { 
			llhoodS = log(lhoodS1) + log(lhoodS2);
		} 

		
	/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */

	llstat = 2.0 * (llhoodS - llhoodP);
	out->width(14);
	out->fill('x');
	/* Significance at the 0.05, 0.01, 0.001 levels requires the statistic, with 1 degrees of freedom, to exceed 3.841, 6.635, and 10.827, respectively. */
	*out << std::fixed  <<std::setprecision(7) << pro.getids() << '\t' << pro.getname(0) << '\t' << pro.getname(1) <<'\t' << pml << '\t' << pmlA << '\t' << pmlB << '\t' << eml << '\t'<< pro.getcoverage(a) << '\t' << pro.getcoverage(b) << '\t' << llstat << std::endl;
	}
exit(0);


	std::cout << t << std::endl;
};
