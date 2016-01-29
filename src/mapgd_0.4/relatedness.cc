#include "relatedness.h"

#define BUFFER_SIZE 500

std::map <Genotype_pair, size_t> hash_genotypes (Indexed_file <population_genotypes> &gcf_in, size_t x, size_t y)
{
	population_genotypes genotypes;
//	gcf_in.?;
	std::map <Genotype_pair, size_t> counts;
	while(gcf_in.is_open() ){
		gcf_in.read(genotypes);
		counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 2)]+=1;
	}
	return counts;
}

/*Does a regression of allele frequency of the samples on the popualtion allele frequency*/
void set_e(Relatedness &relatedness, std::map <Genotype_pair, size_t> &hashed_genotypes)
{
}

/*Guess starting values of relatedness for the maximization procedure*/
void gestimate(Relatedness &relatedness, std::map <Genotype_pair, size_t> &hashed_genotypes)
{
}

/*Maximizes the relatedness*/
void maximize(Relatedness &relatedness, std::map <Genotype_pair, size_t> &hashed_genotypes)
{
}


int estimateRel(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

/*	std::string infile="";
	std::vector <size_t> ind;*/
	std::string gcf_name="", rel_name="";

	env_t env;
	env.setname("mapgd relatedness");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	&gcf_name,	&arg_setstr, 	"an error occured while displaying the help message.", "input file name (default stdout)");
	env.optional_arg('o', "output", &rel_name,	&arg_setstr, 	"an error occured while displaying the help message.", "output file name (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <population_genotypes> gcf_in; 	// Open the file with genotypic probabilities.
	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;			//The class to read to
	population_genotypes genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), std::fstream::in);
	else
		gcf_in.open(std::fstream::in);

	if (rel_name.size()!=0)
		rel_out.open(rel_name.c_str(), std::fstream::out);
	else
		rel_out.open(std::fstream::out);

	genotype=gcf_in.read_header();			//This gives us the sample names.
	rel_out.write_header(relatedness);

	std::map <Genotype_pair, size_t> hashed_genotypes;

	for (size_t x=0; x<relatedness.size(); ++x){
		for (size_t y=x+1; x<relatedness.size(); ++x){
			relatedness.set_X_name(genotype.get_sample_names()[x]);
			relatedness.set_Y_name(genotype.get_sample_names()[y]);
			hashed_genotypes=hash_genotypes(gcf_in, x, y);
			set_e(relatedness, hashed_genotypes);
			gestimate(relatedness, hashed_genotypes);
			maximize(relatedness, hashed_genotypes);
			rel_out.write(relatedness);
		}
	}

	return 0;					//Since everything worked, return 0!.
}
