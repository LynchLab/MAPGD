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

	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "Prints this message");				//DONE
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "Prints the program version");			//DONE

	env.command(	' ',"allele", 	&estimateInd,	 		"an error occured while calling ei", "Estimates allelel frequencies using individual data");		//DONE
	env.command(	' ',"filter",	&filter,	 		"an error occured while calling filter", "Filter sites in '.map' files");			 	//DONE
	env.command(	' ',"genotype",	&map2genotype,	 		"an error occured while calling genotype", "Calculate genotype probabilities for individuals"); 	//DONE
//	env.command(	' ',"vcf",	&vcf,		 		"an error occured while calling vcf", "Convert output to vcf format.");					//TODO
	env.command(	' ',"linkage", 	&PopLD,	 			"an error occured while calling fst", "Estimates linkage disequilibrium between loci");			//TODO
	env.command(	' ',"pool",	&estimatePooled, 		"an error occured while calling ep", "Estimates allele frequencies using pooled data");			//DONE
	env.command(	' ',"proview", 	&proview, 			"an error occured while calling proview", "Prints data in the '.pro' file quartet format");		//DONE
	env.command(	' ',"relatedness",&estimateRel,	 		"an error occured while calling er", "Estimates the pairwise relatedness of individuals");		//TODO
	env.command(	' ',"sam2idx",	&sam2idx,	 		"an error occured while calling genotype", "Reformats a sam header file to a idx file"); 		//DONE

//	env.command(	' ',"read", 	&readsql, 			"an error occured while calling proview", "Reads data from the SQL database");				//TODO
	env.command(	' ',"write", 	&writesql, 			"an error occured while calling proview", "Writes data to the SQL database");				//TODO

//	env.command(	' ',"mlrho", 	&mlrho, 			"an error occured while calling proview", "Runs mlrho");						//TODO
	

//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");
//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");
//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) != 0 ) printUsage(env);
	return 0;
}
