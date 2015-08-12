#include <iostream>
#include <cstdio>
#include <cstring>

#include "interface.h"
#include "commands.h"

using namespace std;

/*! \brief Our main function. Parses commandline arguments etc.
 *	returns 0 iff seccessful
 */
int main (int argc, char* argv[])
{
	env_t env;
	env.setname("mapgd");
	env.setver("1.2");
	env.setauthor("Matthew Ackerman, Bernard Haubold, Michael Lynch, and Takahiro Maruki.");
	env.setdescription("A program for maximum-likelihood analysis of population genomic data.\n Please direction questions to matthew.s.ackerman@gmail.com");

	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");

	env.command(	' ',"write", 	&write_command,			"an error occured while calling write", "write data in a format the humans can use");
	env.command(	' ',"read", 	&read_command,			"an error occured while calling read", "stream data in a format mapgd commands can use");
	env.command(	' ',"ep",	&estimate_pooled, 		"an error occured while calling ep", "estimates allele frequencies using pooled data");
	env.command(	' ',"cp", 	&compare_pooled, 		"an error occured while calling cp", "compares allele frequencies between populations using pooled data");

	env.command(	' ',"ei", 	&estimate_individual, 		"an error occured while calling ei", "estimates allelel frequencies using individual data");

	env.command(	' ',"fst", 	&estimate_fst,	 		"an error occured while calling fst", "estimates fst between two populations");
	
//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");

//	env.command(	' ',"er", 	&estimateRel,	 		"an error occured while calling er", "calculates the relatedness of individuals within a popualtions.");

//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");

//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) != 0 ) printUsage(env);

	return 0;
}
