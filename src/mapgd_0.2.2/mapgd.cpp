#include <iostream>
#include <cstdio>
#include <cstring>

#include "interface.h"
#include "proview.h"
#include "comparePooled.h"
#include "estimatePooled.h"
#include "estimateInd.h"
//#include "calcInd.h"
//#include "estimateRel.h"

using namespace std;

/** @brief Our main function.
  * Parses commandline arguments etc.
  * @returns 0 iff seccessful
**/
int main (int argc, char* argv[])
{
	env_t env;
	env.setname("mapgd");
	env.setver("2.1");
	env.setauthor("Matthew Ackerman, Bernard Haubold, Michael Lynch, and Takahiro Maruki.");
	env.setdescription("A program for maximum-likelihood analysis of population genomic data.\n Please direction questions to matthew.s.ackerman@gmail.com");

	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");

	env.command(	' ',"proview", 	&proview, 			"an error occured while calling proview", "prints data in the '.pro' file quartet format");
	env.command(	' ',"ep", 	&estimatePooled, 		"an error occured while calling ep", "estimates allele frequencies using pooled data");
	env.command(	' ',"cp", 	&comparePooled, 		"an error occured while calling cp", "compares allele frequencies between populations using pooled data");

	env.command(	' ',"ei", 	&estimateInd,	 		"an error occured while calling ei", "estimates allelel frequencies using individual data");
	
//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");

//	env.command(	' ',"er", 	&estimateRel,	 		"an error occured while calling er", "calculates the relatedness of individuals within a popualtions.");

//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");

//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) ) printUsage(env);

	printUsage(env);

	exit(0);
}
