#include "mpi_relatedness.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

#define MASTER 0 

int estimateRel(int argc, char *argv[])
{
	std::cerr << "Warning: this program should generate AICc's. However, it doesn't and its _ll variables don't mean that much.\n"; 

	/* All the variables that can be set from the command line */

	std::string gcf_name="", rel_name="";

	Environment env;
	env.set_name("mapgd relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	gcf_name, "an error occurred while displaying the help message.", "input filename (default stdout)");
	env.optional_arg('o', "output", rel_name, "an error occurred while displaying the help message.", "output filename (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;		//The class to read to
	Population genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), std::ios::in);
	else
		gcf_in.open(std::ios::in);

	if (rel_name.size()!=0)
		rel_out.open(rel_name.c_str(), std::ios::out);
	else
		rel_out.open(std::ios::out);
	
	genotype=gcf_in.read_header();			//This gives us the sample names.

	gcf_mem.open(&file_buffer, std::ios::out );
	gcf_mem.set_index(gcf_in.get_index() );
	gcf_mem.write_header(genotype);

	while(gcf_in.table_is_open() ){
		gcf_in.read(genotype);
		gcf_mem.write(genotype);
	}
	gcf_mem.close();

	rel_out.write_header(relatedness);

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	std::map <Genotype_pair_tuple, size_t> down_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();
	
	Relatedness buffer_rel[BUFFER_SIZE];        
	std::fill_n(buffer_rel, BUFFER_SIZE, relatedness);

	int taskid,	/* task ID - also used as seed number */
	numtasks,	/* number of tasks */
	rc,             /* return code */
	i;		/* new comment */

	MPI_Status status;
	/* Obtain number of tasks and task ID */
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

#if DEBUG
	printf ("MPI task %d has started...\n", taskid);
#endif

	size_t done=0;
	for (size_t x=0; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
			z=x*sample_size+(y-x-1)-done;
			if (z>BUFFER_SIZE){
				if (taskid == MASTER) 
				{
					MPI_Waitall(numtasks-1, request, status);
					rel_out.write(buffer_rel[z]);
				}
				done+=BUFFER_SIZE;
				z-=BUFFER_SIZE;
			} 
			if (z%numtasks==taskid){
				relatedness.set_X_name(genotype.get_sample_names()[x]);
				relatedness.set_Y_name(genotype.get_sample_names()[y]);
				hashed_genotypes=hash_genotypes(file_buffer, x, y);
				down_genotypes=downsample_genotypes(file_buffer, x, y);
				relatedness.zero();
				set_e(relatedness, hashed_genotypes);
			//	gestimate(relatedness, hashed_genotypes);
#ifdef EIGEN
				newton(relatedness, down_genotypes);
#else
				maximize(relatedness, down_genotypes);
#endif
				get_llr(relatedness, hashed_genotypes);
				buffer_rel[z]=relatedness;
			}
		}
	}
	rel_out.close();
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

