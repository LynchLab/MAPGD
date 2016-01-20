#include <iostream>
#include <cstdio>
#include <cstring>

#include "interface.h"
#include "proview.h"
#include "commands.h"


using namespace std;

/*! \brief Our main function. Parses commandline arguments etc.
 *	returns 0 iff seccessful
 */
int main (int argc, char* argv[])
{
#if defined(DEBUG)
	std::cerr << "Warning: This version of mapgd has been compiled in debug mode, and will run slowly.\n";
#endif
	env_t env;
	env.setname("mapgd");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman, Bernard Haubold, Michael Lynch, and Takahiro Maruki.");
	env.setdescription("A program for maximum-likelihood analysis of population genomic data.\n Please direct questions to matthew.s.ackerman@gmail.com");

	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "Prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "Prints the program version");

	env.command(	' ',"cp", 	&comparePooled, 		"an error occured while calling cp", "Compares allele frequencies between populations using pooled data");

	env.command(	' ',"allele", 	&estimateInd,	 		"an error occured while calling ei", "Estimates allelel frequencies using individual data");
	env.command(	' ',"convert", 	&convert,	 		"an error occured while calling convert", "Converts plain text genotype files into binary format.");

//	env.command(	' ',"genotype",	&genotype,	 		"an error occured while calling genotype", "Calculate genotype probabilities for individuals.");
//	env.command(	' ',"filter",	&filter,	 		"an error occured while calling filter", "Estimates linkage disequilibrium between loci.");
//	env.command(	' ',"vcf",	&vcf,		 		"an error occured while calling vcf", "Convert output to vcf format.");

	env.command(	' ',"linkage", 	&PopLD,	 			"an error occured while calling fst", "Estimates linkage disequilibrium between loci");
	env.command(	' ',"pool",	&estimatePooled, 		"an error occured while calling ep", "Estimates allele frequencies using pooled data");
	env.command(	' ',"proview", 	&proview, 			"an error occured while calling proview", "Prints data in the '.pro' file quartet format");

//	env.command(	' ',"relatedness",&estimateRel,	 		"an error occured while calling er", "Estimates the pairwise relatedness of individuals");
//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");

//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");

//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) != 0 ) printUsage(env);
	return 0;
}
