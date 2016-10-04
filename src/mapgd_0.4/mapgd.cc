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
#if defined(DEBUG)
	std::cerr << "Warning: This version of mapgd has been compiled in debug mode, and will run slowly.\n";
#endif
	Environment env;
	env.set_name("mapgd");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman, Bernard Haubold, Michael Lynch, and Takahiro Maruki.");
	env.set_description("A program for maximum-likelihood analysis of population genomic data. Please direct questions to matthew.s.ackerman@gmail.com");

	env.flag(	'a',"avail", 	&env, 		&flag_commands, "an error occurred while displaying the help message", "Prints a list of available commands");		//DONE
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "Prints this message");				//DONE
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "Prints the program version");		//DONE
	env.command(	' ',"allele", 	&estimateInd,	 		"an error occurred while calling allele", "Estimates allele frequencies using individual data");	//DONE
	env.command(	' ',"filter",	&filter,	 		"an error occurred while calling filter", "Filter sites in '.map' files");			 	//DONE
	env.command(	' ',"filterpool",&filter_pool,	 		"an error occurred while calling filter", "Filter sites in '.map' files");			 	//DONE
	env.command(	' ',"genotype",	&map2genotype,	 		"an error occurred while calling genotype", "Calculate genotype probabilities for individuals"); 	//DONE
//	env.command(	' ',"vcf",	&vcf,		 		"an error occured while calling vcf", "Convert output to vcf format.");					//TODO
	env.command(	' ',"linkage", 	&PopLD,	 			"an error occurred while calling linkage", "Estimates linkage disequilibrium between loci");		//DONE
	env.command(	' ',"pool",	&estimatePooled, 		"an error occurred while calling pool", "Estimates allele frequencies using pooled data");		//DONE
	env.command(	' ',"proview", 	&proview, 			"an error occurred while calling proview", "Prints data in the '.pro' file quartet format");		//DONE
	env.command(	' ',"relatedness",&estimateRel,	 		"an error occurred while calling relatedness", "Estimates the pairwise relatedness of individuals");	//DONE
	env.command(	' ',"sam2idx",	&sam2idx,	 		"an error occurred while calling sam2idx", "Reformats a sam header file to a idx file"); 		//DONE
	env.command(	' ',"vcf",	&make_vcf,	 		"an error occurred while calling vcf", "Prints as a vcf file"); 		//DONE
	env.command(	' ',"help", 	&mapgd_help, 			"an error occurred while calling help", "Prints helpful information");					//DONE

#ifndef NOSQL
	env.command(	' ',"read", 	&readsql, 			"an error occurred while calling read", "Reads data from the SQL database");				//TODO
	env.command(	' ',"write", 	&writesql, 			"an error occurred while calling write", "Writes data to the SQL database");				//DONE
#endif
	env.set_footer(" 'mapgd help [command]' will print the manual page for that command, e.g. mapdg help genotype will show general information the genotype command. mapgd help [COLUMN NAME] will \
give information about the column names that appear in the first two lines of files.");
//	env.command(	' ',"mlrho", 	&mlrho, 			"an error occured while calling proview", "Runs mlrho");						//TODO
//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");
//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");
//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);
	return 0;
}
