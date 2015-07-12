/* 

Program estimatePooled.cpp:

	1) input from a list of site-specific quartets from a pooled population sample;

	2) identify the major and minor alleles, obtain the maximum-likelihood estimates of allele frequencies, and significance levels.

The designated major nucleotide is simply the one with the highest rank, and the minor nucleotide the one with the second highest rank.
	If the top three ranks are all equal, the site is treated as unresolvable, with both major and minor nucelotides designated in the output by a *.
	If the second and third ranks are equal but lower than the major-nucleotide count, the site is treated as monomorphic, with the minor nucleotide designated by a *.
	
Input File: six columns; first two are arbitrary identifiers of the site (e.g., chromosome and position); final four are integer numbers of reads of A, C, G, T.
	Columns are tab delimited.
	Default name is "datain.txt".

Output File: two columns of site identifiers; major allele; minor allele; major-allele frequency; minor-allele frequency; error rate; likeilhood-ratio test of polymorphism. 
	Columns are comma delimited.
	Default name is "dataout.txt".

***** Note that significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

***** The 95% support interval can be obtained by determining the changes in the estimate of p in both directions required to reduce the log likelihood by 3.841
      although this is not currently implemented. 
*/

#include "estimate-pooled.h"

int estimatePooled(int argc, char *argv[])
{

	/* variables that can be set from the command line */
	std::string infile="";
	std::string outfile="";

	bool verbose=false;
	bool quite=false;
	float_t EMLMIN=0.00;
	float_t a=0.00;

	std::vector <int> pop;

	/* sets up the help messages and options */

	env_t env;
	env.setname("mapgd ep");
	env.setver("2.0");
	env.setauthor("Michael Lynch and Matthew Ackerman");
	env.setdescription("Uses a maximum likelihood approach to estimates allele frequencies in pooled population genomic data");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file", "sets the input file for the program (default 'datain.txt')");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file", "sets the output file for the program (default 'dataout.txt')");
	env.optional_arg('p',"populations", &pop, &arg_setvectorint, "please provide a list of integers", "choose populations to compare");
	env.optional_arg('m',"minerror", &EMLMIN, &arg_setfloat_t, "please provide a float", "minimum error rate.");
	env.optional_arg('a',"alpha", &a, &arg_setfloat_t, "please provide a float", "alpha.");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'V',"verbose", 	&verbose,	&flag_set, 	"an error occured", "prints more information while the command is running");
	env.flag(	'q',"quite", 	&quite,		&flag_set, 	"an error occured", "prints less information while the command is running");

	/* gets any command line options */
	if ( parsargs(argc, argv, env) ) printUsage(env);
	/* Point to the input and output files. */


	std::ostream *out;
	std::ofstream outFile;
	out=&std::cout;
	profile pro;

	if (infile.size()!=0) {if (pro.open(infile.c_str(), "r")==NULL) {printUsage(env);} }
	else pro.open("r");

	float_t pmaj, maxll;						/* fraction of putative major reads among the total of major and minor */
	float_t popN, eml, enull, llhoodMS, llhoodPS, llhoodFS;		/* maximum-likelihood estimates of the error rate, null error rate, and major-allele frequency */
	float_t *llhoodM, *llhoodP, *llhoodF, *llstat, *llstatc,*pmlP;	/* loglikelihood sums */

	/* Open the output and input files. */

	if (outfile.size()!=0) {
		outFile.open(outfile, std::ofstream::out);
		if (!outFile) printUsage(env);
		out=&outFile;
	};


	pro.maskall();

	if ( pop.size()==0 ) { 
		pop.clear();
		for (count_t x=0; x<pro.size(); ++x) pop.push_back(x);
	};

	for (count_t x=0; x<pop.size(); ++x) pro.unmask(pop[x]);

	*out << "id1\t\tid2\tmajor\tminor\t";
	for (int x=0; x<pop.size(); ++x) *out << "Freq_P\tll_poly\tll_fixed\tcov\t";
	*out << "\tTotcov\tError" << std::endl;
	/* Start inputting and analyzing the data line by line. */
	/* The profile format can contain quartets for many populations and these populations need to be
	 * iterated across.
	 */

	llhoodP=new float_t[pop.size()];
	llhoodF=new float_t[pop.size()];
	llhoodM=new float_t[pop.size()];
	pmlP=new float_t[pop.size()];
	llstat=new float_t[pop.size()];
	llstatc=new float_t[pop.size()];	//

	lnmultinomial multi(4);			//a class for our probability function (which is just a multinomial distribution).
	allele_stat site;
	site_t line;

	while (pro.read(line)!=EOF ){
		//sorts the reads at the site from most frequenct (0) to least frequenct (3) (at metapopulation level);

		/* Calculate the ML estimates of the major / minor allele frequencies and the error rate. */
                line.sort();

		site.major=line.getindex(0);
		site.minor=line.getindex(1);

                /* 3) CALCULATE THE LIKELIHOOD OF THE POOLED POPULATION DATA. */

                /* Calculate the ML estimates of the major pooled allele frequency and the error rate. */
		popN=line.getcoverage();

                site.error= (float_t) ( popN - line.getcount(0) ) / ( (float_t) popN );
		site.error= (float_t) (line.getcoverage()-line.getcount(0)-line.getcount(1) )*3. / ( 2.*(float_t) popN ) ;

		if (site.error<EMLMIN) site.error=EMLMIN;

		multi.set(&monomorphicmodel, site);

                for (int x=0; x<pop.size(); ++x) llhoodM[x]=multi.lnprob(line.getquartet(pop[x]) );

		site.error= (float_t) (line.getcoverage()-line.getcount(0)-line.getcount(1) )*3. / ( 2.*(float_t) popN ) ;

		if (site.error<EMLMIN) site.error=EMLMIN;

                /* Calculate the likelihoods under the reduced model assuming no variation between populations. */

		multi.set(&fixedmorphicmodel, site);
                for (int x=0; x<pop.size(); ++x) llhoodF[x]=multi.lnprob(line.getquartet(pop[x]) );

                /* 4) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF monomorphism */

                for (int x=0; x<pop.size(); ++x) {
                        pmaj = (float_t) line.getcount(pop[x],0) / (float_t) ( line.getcount(pop[x], 1) + line.getcount(pop[x],0) );
			
                        pmlP[x] = (pmaj * (1.0 - (2.0 * eml / 3.0) ) ) - (eml / 3.0);
                        pmlP[x] = pmlP[x] / (1.0 - (4.0 * eml / 3.0) );
			if (pmlP[x]<0) pmlP[x]=0.;
			if (pmlP[x]>1) pmlP[x]=1.;

			site.freq=pmlP[x];
			multi.set(&polymorphicmodel, site);
                        llhoodP[x]=multi.lnprob(line.getquartet(pop[x]) );
                };

                /* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */

                maxll=0;
                for (int x=0; x<pop.size(); ++x){
                        llhoodMS=0;
                        llhoodPS=0;
                        llhoodFS=0;
			/*
                        for (int y=0; y<pop.size(); ++y){
                                llhoodMS+=llhoodM[y];
                                if (x!=y){
					llhoodPS+=llhoodM[y];
					llhoodFS+=llhoodM[y];
				}
                                else {
					llhoodPS+=llhoodP[y];
					llhoodFS+=llhoodF[y];
				};
                        };*/
                        llstat[x] = (2.0 * (llhoodP[x] - llhoodM[x]) );
                        llstatc[x] = (2.0 * (llhoodP[x] - llhoodF[x]) );
                        if (llstat[x]>maxll) maxll=llstat[x];
                };

                if (maxll>=a){
                        *out << std::fixed << std::setprecision(7) << pro.getids(line) << '\t' << line.getname(0) << '\t' << line.getname_gt(1) << '\t';
                        for (int x=0; x<pop.size(); ++x){
                                if ( line.getcoverage(pop[x])==0) *out << std::fixed << std::setprecision(7) << "NA" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\t';
                                else *out << std::fixed << std::setprecision(7) << pmlP[x] <<'\t' << llstat[x] << '\t' << llstatc[x]<< '\t' << line.getcoverage(pop[x]) <<'\t';
                        };
                        /* Significance at the 0.05, 0.01, 0.001 levels requires the statistic, with 1 degrees of freedom, to exceed 3.841, 6.635, and 10.827, respectively. */
                        *out << line.getcoverage() << '\t' << site.error << std::endl;
                };

	}
        delete [] llhoodP;
        delete [] llhoodF;
        delete [] llhoodM;
        delete [] pmlP;
        delete [] llstat;
        delete [] llstatc;
	if (outFile.is_open()) outFile.close();
	/* Ends the computations when the end of file is reached. */
	exit(0);
}
