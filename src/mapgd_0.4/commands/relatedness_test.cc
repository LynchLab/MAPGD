#include "relatedness_test.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

/*Guess starting values of relatedness for the maximization procedure*/

int testRel(int argc, char *argv[])
{
	/* All the variables that can be set from the command line */

	bool l2o=false;
	std::string gcf_name="", rel_name1="", rel_name2="";

	Environment env;
	env.set_name("mapgd trelatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Test for significant differences between relatedness estimates.");

	env.positional_arg('i', "input",gcf_name, "an error occurred while displaying the help message.", "input filename (default stdin)");
	env.positional_arg('1', "rel1", rel_name1, "an error occurred while displaying the help message.", "input filename (default stdin)");
	env.positional_arg('2', "rel2", rel_name2, "an error occurred while displaying the help message.", "input filename (default stdin)");
	env.flag(	'l', "l2o", 	&l2o, 		&flag_set, 	"an error occurred while displaying the help message.", "attempt to use leave two out (not working");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_in1; 		// Open a file for input.
	Flat_file <Relatedness> rel_in2; 		// Open a file for input.

	Relatedness rel2;			//The class to read to
	Relatedness rel1;			//The class to read to
	Population genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), std::ios::in);
	else
		gcf_in.open(std::ios::in);

	if (rel_name1.size()!=0)
		rel_in1.open(rel_name1.c_str(), std::ios::in);
	else
		rel_in1.open(std::ios::in);

	if (rel_name2.size()!=0)
		rel_in2.open(rel_name2.c_str(), std::ios::in);
	else
		rel_in2.open(std::ios::in);
	
	genotype=gcf_in.read_header();			//This gives us the sample names.

	gcf_mem.open(&file_buffer,  BINARY | WRITE  );
	gcf_mem.set_index(gcf_in.get_index() );
	gcf_mem.write_header(genotype);

	while(gcf_in.table_is_open() ){
		gcf_in.read(genotype);
		gcf_mem.write(genotype);
	}
	gcf_mem.close();

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();

	float_t r1, r2;

	rel1=rel_in1.read_header();			//This gives us the sample names.
	rel2=rel_in2.read_header();			//This gives us the sample names.

	while(rel_in1.table_is_open() ){
		rel_in1.read(rel1);
		rel_in2.read(rel2);

		if (rel1.X_!=rel2.X_ || rel1.Y_!=rel2.Y_){
			std::cerr << __FILE__ << ":" << __LINE__ << "relatedness in files do not compare the same individual. exiting.\n";
			break;
		}

		hashed_genotypes=downsample_genotypes(file_buffer, rel1.X_, rel1.Y_, l2o);
		std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );

		gsl_vector *x;
		
		x=gsl_vector_alloc(7);

		gsl_vector_set(x, 0, rel1.f_X_);
		gsl_vector_set(x, 1, rel1.f_Y_);
		gsl_vector_set(x, 2, rel1.theta_XY_);
		gsl_vector_set(x, 3, rel1.gamma_XY_);
		gsl_vector_set(x, 4, rel1.gamma_YX_);
		gsl_vector_set(x, 5, rel1.Delta_XY_);
		gsl_vector_set(x, 6, rel1.delta_XY_);

		r1=rel_ll(x, &hashed_genotypes_vector);	

		gsl_vector_set(x, 0, rel2.f_X_);
		gsl_vector_set(x, 1, rel2.f_Y_);
		gsl_vector_set(x, 2, rel2.theta_XY_);
		gsl_vector_set(x, 3, rel2.gamma_XY_);
		gsl_vector_set(x, 4, rel2.gamma_YX_);
		gsl_vector_set(x, 5, rel2.Delta_XY_);
		gsl_vector_set(x, 6, rel2.delta_XY_);

		r2=rel_ll(x, &hashed_genotypes_vector);
		if (r1>r2) std::cout << "r1 < r2 : " << r1-r2 << ", " << -r2 << std::endl;
		else std::cout << "r2 < r1 : " << r2-r1 << ", " << -r1 << std::endl;
	}
	return 0;					//Since everything worked, return 0!.
}

#else 

int 
testRel(int argc, char *argv[])
{
	std::cerr << "This command depends on gsl, which could not be found. Please e-mail matthew.s.ackerman@gmail.com for help.\n";
	return 0;
}

#endif 

