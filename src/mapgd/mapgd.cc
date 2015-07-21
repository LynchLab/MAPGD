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
	env_t env;
	env.setname("mapgd");
	env.setver("1.2");
	env.setauthor("Matthew Ackerman, Bernhard Haubold, Michael Lynch, and Takahiro Maruki.");
	env.setdescription("A program for maximum-likelihood analysis of population genomic data.\n Please direction questions to matthew.s.ackerman@gmail.com");

	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message.");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version.");

	env.command(	' ',"proview", 	&proview, 			"an error occured while calling proview", "prints data in the '.pro' file quartet format.");

//#	env.command(	' ',"gcfview", 	&gcfview, 			"an error occured while calling gcfview", "prints data in the '.gcf' file format.");

	env.command(	' ',"ep", 	&estimatePooled, 		"an error occured while calling ep", "estimates allele frequencies using pooled data.");
	env.command(	' ',"cp", 	&comparePooled, 		"an error occured while calling cp", "compares allele frequencies between populations using pooled data.");

	env.command(	' ',"ei", 	&estimateInd,	 		"an error occured while calling ei", "estimates allele frequencies using individual data.");
	
	if ( parsargs(argc, argv, env) ) printUsage(env);
	return 0;
}
