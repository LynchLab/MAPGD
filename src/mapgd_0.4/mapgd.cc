#include <iostream>
#include <cstdio>
#include <cstring>

#include "interface.h"
#include "commands.h"


using namespace std;

/*! \brief Our main function. Parses command line arguments etc.
 *	returns 0 iff successful
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
	env.flag(	'O',"opts", 	&env, 		&flag_options, "an error occurred while displaying the help message", "Prints a list of available options");		//DONE
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "Prints this message");				//DONE
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "Prints the program version");		//DONE
	env.command(	' ',"allele", 	&allele_cmd,	 		"an error occurred while calling allele", "Estimates allele frequencies using individual data");	//DONE
	env.command(	' ',"fastview",	&fastview,	 		"an error occurred while calling filter", "Quickly displays contents of a file");			 	//DONE
//	env.command(	' ',"filter",	&meta_filter,	 		"an error occurred while calling filter", "Filter contents of a file ");			 	//DONE
	env.command(	' ',"filter",	&filter,	 		"an error occurred while calling filter", "Filter sites in '.map' files");			     	//DONE
	env.command(	' ',"filterpool",&filter_pool,	 		"an error occurred while calling filter", "Filter sites in '.pol' files");			 	//DONE
//	env.command(	' ',"filterpro", &filter_pro,	 		"an error occurred while calling filter", "Filter sites in '.pro' files");			 	//DONE
	env.command(	' ',"filtergcf",&filter_genotype,	 	"an error occurred while calling filter", "Filter sites in '.gcf' files");			 	//DONE
	env.command(	' ',"genotype",	&map2genotype,	 		"an error occurred while calling genotype", "Calculate genotype probabilities for individuals"); 	//DONE
	env.command(	' ',"linkage", 	&linkage_disequilibrium,	 			"an error occurred while calling linkage", "Estimates linkage disequilibrium between loci");		//DONE
	env.command(	' ',"pool",	&estimatePooled, 		"an error occurred while calling pool", "Estimates allele frequencies using pooled data");		//DONE
	env.command(	' ',"proview", 	&proview, 			"an error occurred while calling proview", "Prints data in the '.pro' file quartet format");		//DONE
	env.command(	' ',"relatedness",&estimateRel,	 		"an error occurred while calling relatedness", "Estimates the pairwise relatedness of individuals");	//DONE
	env.command(	' ',"reltest", &testRel,	 		"an error occurred while calling relatedness", "Test for sig dif between relatedness estimates");	//DONE
	env.command(	' ',"sam2idx",	&sam2idx,	 		"an error occurred while calling sam2idx", "Reformats a sam header file to a idx file"); 		//DONE
	env.command(	' ',"keyinfo",	&test_keys,	 		"an error occurred while calling sam2idx", "Displays information regarding keys (i.e. column names)"); 		//DONE
//	env.command(	' ',"simulate",	&simulate,	 		"an error occurred while calling sam2idx", "Test, do not use"); 		//DONE

#ifndef NOHTS
	env.command(	' ',"writevcf",	&make_vcf,	 		"an error occurred while calling vcf", "Prints as a vcf file"); 		//DONE
	env.command(	' ',"writevcf2",&make_vcf2,	 		"an error occurred while calling vcf", "Prints as a vcf file"); 		//DONE
	env.command(	' ',"readvcf",	&read_vcf,	 		"an error occurred while calling vcf", "Reads a vcf file"); 		//DONE
#endif
	env.command(	' ',"readbed",	&read_bed,	 		"an error occurred while calling vcf", "Reads a bed file"); 		//DONE
	env.command(	' ',"readphen",	&read_pheno,	 		"an error occurred while calling vcf", "Reads plink's pheno file"); 		//DONE
	env.command(	' ',"help", 	&mapgd_help, 			"an error occurred while calling help", "Prints helpful information");					//DONE

#ifndef NOSQL
	env.command(	' ',"read", 	&readsql, 			"an error occurred while calling read", "Reads data from the SQL database");				//TODO
	env.command(	' ',"write", 	&writesql, 			"an error occurred while calling write", "Writes data to the SQL database");				//DONE
#endif
	env.set_footer(" '\x1b[1mmapgd help\x1b[0m \x1b[4mCOMMAND\x1b[0m' will print the manual page for that command, e.g. \x1b[1mmapdg help genotype\x1b[0m will show general information for the genotype command. \x1b[1mmapgd help\x1b[0m \x1b[4mCOLUMN NAME\x1b[0m will \
give information about the column names that appear in the first two lines of files.");
//	env.command(	' ',"mlrho", 	&mlrho, 			"an error occured while calling proview", "Runs mlrho");						//TODO
//	env.command(	' ',"calc", 	&calcInd,	 		"an error occured while calling calc", "print log likelihood of ...");
//	env.command(	' ',"ci", 	&compareIndividual, 		"an error occured while calling ci", "compares allele frequencines between population using indivdual data ");
//	env.command(	' ',"eq", 	&estimateQuality, 		"an error occured while calling eq", "estimates error rate as a function of illumina quality score");

	if ( parsargs(argc, argv, env) != 0 ) print_usage(env);
	return 0;
}
