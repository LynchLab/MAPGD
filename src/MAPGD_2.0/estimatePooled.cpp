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

#include "estimatePooled.h"
#define EMLMIN	0.00001

int estimatePooled(int argc, char *argv[])
{

/* variables that can be set from the command line */
	std::string infile="datain.txt";
	std::string outfile="dataout.txt";

	bool verbose=false;
	bool quite=false;
	int p=1;

/* sets up the help messages and options */

	env_t env;
	env.setname("mapgd ep");
	env.setver("1.0");
	env.setauthor("Michael Lynch");
	env.setdescription("Uses a maximum likelihood approach to estimates allele frequencies in pooled population genomic data");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file", "sets the input file for the program (default 'datain.txt')");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file", "sets the output file for the program (default 'dataout.txt')");
	env.optional_arg('p',"output", 	&p,		&arg_setint, 	"an error occured while setting p", "the population to analize");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'V',"verbose", 	&verbose,	&flag_set, 	"an error occured", "prints more information while the command is running");
	env.flag(	'q',"quite", 	&quite,		&flag_set, 	"an error occured", "prints less information while the command is running");

	/* gets any command line options */
	if ( parsargs(argc, argv, env) ) printUsage(env);
	p--;
	/* Point to the input and output files. */

	FILE *outstream=NULL;
	std::ostream *out;
	out=&std::cout;
	profile pro;

	if (infile.size()!=0) {if (pro.open(infile.c_str(), 'r')==NULL) {printUsage(env);} }
	else pro.open('r');

	int coverage;						/* total number of reads at the site */

	int nmax;						/* number of reads for the most abundant type (putative major allele) */
	int nmin;						/* number of reads for the second most abundant type (putative minor allele) */
	int ner1;						/* number of reads for the third highest ranked nucleotide */

	int major, minor;					/* identities of major and minor alleles */

	double pmaj;						/* fraction of putative major reads among the total of major and minor */
	double eml, pml;					/* maximum-likelihood estimates of the error rate and major-allele frequency */

	double lphimaj, lphimin, lphie;				/* terms for the likelihood expression */
	double emltheta;
	double llhoodp, llhoodm;				/* log likelihoods under the assumption of polymorphism and monomorphism */

	int nerrors, nerrorsm;					/* errors inferred in the full- and reduced-model likelihood analyses */

	double errorm;						/* estimated error rate in the reduced model */
	double llstat;						/* difference in the log likelihoods */

	bool fliped;

	/* Open the output and input files. */

	outstream = fopen(outfile.c_str(), "w");

	if (outstream==NULL) {std::cerr << "failed to open file " << outfile << " for writing\n"; printUsage(env);}

	*out << "id1\t\tid2\tmajor\tminor\tpml\teml\tcov\tllstat" << std::endl;
	/* Start inputting and analyzing the data line by line. */

	while (pro.read()!=EOF ){
		/* The profile format can contain quartets for many populations and these populations need to be
		 * iterated accors
		 */
		//sorts the reads at the site from most frequenct (0) to least frequenct (3) (at metapopulation level);
		pro.maskall();
		pro.unmask(p);
		pro.sort();
		major=0;
		minor=1;
		nmax=pro.getcount(0);
		//set the minor allele count to the second most frequent read (in the meta population);
		nmin=pro.getcount(1);
		//set the first error  to the third most frequent read (in the meta population);
		ner1=pro.getcount(2);

		if  ( (nmax == nmin) && (nmin == ner1) ) {
			pml = 0.0;
			eml = 0.0;
			llstat = 0.0; 
			major = 4;
			minor = 4;
		}

		if  ( (nmax > nmin) && (nmin == ner1) ) {
			pml = 1.0;
			llstat = 0.0; 
			minor = 4;
		}
		
		*out << pro.getids() << '\t' << pro.getname(major) << '\t' << pro.getname(minor) << '\t'; 

		for (int s=0; s<pro.samplesize(); ++s){
			fliped=false;
			//get the coverage for the sth population.
			coverage = pro.getcoverage(s);
			/* Identify the counts for the putative major and minor alleles.  */
			
			//set the major allele count to the most frequent read (in the meta);
			nmax=pro.getcount(s,0);
			//set the minor allele count to the second most frequent read (int the meta);
			nmin=pro.getcount(s,1);

			if (nmax<nmin) {std::swap(nmax, nmin); fliped=true;}
			/* Calculate the ML estimates of the major / minor allele frequencies and the error rate. */

			pmaj = ( (double) nmax) / ( ((double) nmax) + ((double) nmin) );

			eml = 1.5 * ( ((double) coverage) - ((double) nmax) - ((double) nmin) ) / ((double) coverage);

			//if(eml==0) eml=EMLMIN;
	
			pml = (pmaj * (1.0 - (2.0 * eml / 3.0)) ) - (eml / 3.0);
			pml = pml / (1.0 - (4.0 * eml / 3.0) );

			//if(pml>1.0) pml=0.99999;
			/* Calculate the log likelihoods under: 1) the full model; and 2) the reduced model assuming no polymorphism (major-allele frequency = 1.0). */

			/* First, get the log likelihood under the full model, which assumes a polymorphism. */
	
			emltheta = 1.0 - (4.0 * eml / 3.0);
			lphimaj = (pml * emltheta) + (eml / 3.0);
			lphimin = 1.0 - eml - (pml * emltheta);
			lphie =  eml / 3.0;
		
			nerrors = coverage - nmax - nmin;

			if ( (lphie == 0.0) && (lphimin > 0.0) ) {
				llhoodp = (((double) nmax) * log(lphimaj)) + (((double) nmin) * log(lphimin));}
			else if ( (lphie == 0.0) && (lphimin == 0.0) ) {
				llhoodp = 0.0;}
			else {
				llhoodp = (((double) nmax) * log(lphimaj)) + (((double) nmin) * log(lphimin)) + (((double) nerrors) * log(lphie));}
	

			/* Second, get the log likelihood under the assumption of just one true allele. */

			nerrorsm = coverage - nmax;
		
			errorm = ((double) nerrorsm) / ((double) coverage);
		
			if ( errorm == 0.0 ) {
				llhoodm = 0.0;
			} else 
				llhoodm = ( ((double) nmax) * log(1.0 - errorm) ) + (((double) nerrorsm) * log(errorm/3.0) );
	

			/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */

			llstat = 2.0 * (llhoodp - llhoodm);

			/* If the top three ranks are all equal, flag as unresolvable by denoting the alleles with asterisks. */
	
			/* If the second and third ranks are equal but less than the top, treat as fixed for the putative major allele. */
	
			if (coverage==0){
				eml=0.0;
				llstat=0.0;
			}
			if (!fliped)
				*out <<  pml  << '\t' << eml  << '\t' <<  coverage << '\t' <<  llstat << '\t'; 
			else
				*out <<  (1.0-pml) << '\t' << eml  << '\t' <<  coverage << '\t' <<  llstat << '\t'; 

		}
		*out << std::endl;
	}
	/* Ends the computations when the end of file is reached. */
exit(0);
}
