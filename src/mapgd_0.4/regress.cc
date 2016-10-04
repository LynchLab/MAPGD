#include "relatedness.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

int regress(int argc, char *argv[])
{
	/* All the variables that can be set from the command line */

	std::string ped_name="", rel_name="";

	Environment env;
	env.set_name("mapgd regress");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Performs a least squares regression to obtain population variances.");

	env.optional_arg('r', "input", 	rel_name, "an error occurred while displaying the help message.", "input relatedness (default stdout)");
	env.optional_arg('p', "ped", 	ped_name, "an error occurred while displaying the help message.", "input ped");
	env.optional_arg('o', "output", lsq_name, "an error occurred while displaying the help message.", "output filename (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); 	//Gets all the command line options, and prints usage on failure.

	Flat_file <Relatedness> rel_in; 			// Open a file for output.

	Relatedness relatedness;				//The class to read to

	if (rel_name.size()!=0)
		rel_in.open(rel_name.c_str(), std::ios::in);
	else
		rel_in.open(std::ios::in);

	if (ped_name.size()!=0)
		ped_in.open(ped_name.c_str(), std::ios::out);
	else
		ped_in.open(std::ios::out);
	
	rel_in.read_header(relatedness);

	size_t sample_size=genotype.get_sample_names().size();

	for (size_t x=0; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
			relatedness.set_X_name(genotype.get_sample_names()[x]);
			relatedness.set_Y_name(genotype.get_sample_names()[y]);
		}
	}
	return 0;					//Since everything worked, return 0!.
}

#else 

int 
estimateRel(int argc, char *argv[])
{
	std::cerr << "This command depends on gsl, which could not be found. Please e-mail matthew.s.ackerman@gmail.com for help.\n";
	return 0;
}

#endif 

