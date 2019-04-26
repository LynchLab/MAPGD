#include "surface.h"

#ifdef NOMPI

int surface(int argc, char *argv[])
{
	std::cerr << "Warning: this program should generate AICc's. However, it doesn't and its _ll variables don't mean that much. Also it is way too slow.\n"; 

	/* All the variables that can be set from the command line */

	std::string gcf_name="", rel_name="";

	bool l2o=false, called=false, newton=false;
    int startx=0, starty=1, axis=0;

    std::vector <double> around;
	Environment env;
	env.set_name("mapgd surface");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Prints the likelihood surface for a pair of individuals.");

	env.optional_arg('i', "input", 	gcf_name, "an error occurred while displaying the help message.", "input filename (default stdout)");
	env.optional_arg('x', "x", startx, "an error occurred while displaying the help message.", "individual x (default 0)");
	env.optional_arg('y', "y", starty, "an error occurred while displaying the help message.", "individual y (default 1)");
	env.optional_arg('d', "direction", axis, "an error occurred while displaying the help message.", "individual y (default 1)");
	env.optional_arg(	'a', "around", 	around, "an error occurred while displaying the help message.", "look at points around the point ....");
	env.flag(	'n', "newton", 	&newton, 		&flag_set, "an error occurred while displaying the help message.", "uses the NR optimization. Doesn't work.");
	env.flag(	'l', "l2o", 	&l2o, 		&flag_set, "an error occurred while displaying the help message.", "uses the 'leave 2 out' procedure of calculating allele freq.");
	env.flag(	'c', "called", 	&called,	&flag_set, "an error occurred while displaying the help message.", "ignore genotypic likelihoods and treat data as error free.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Relatedness relatedness;		//The class to read to
	Population genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), READ );
	else
		gcf_in.open( READ );

	
	genotype=gcf_in.read_header();			//This gives us the sample names.

	gcf_mem.open(&file_buffer, BINARY | WRITE );
	gcf_mem.set_index(gcf_in.get_index() );
	gcf_mem.write_header(genotype);

	while(gcf_in.table_is_open() ){
		gcf_in.read(genotype);
		gcf_mem.write(genotype);
	}
	gcf_mem.close();

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	std::map <Genotype_pair_tuple, size_t> down_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();

	relatedness.set_X_name(startx);
	relatedness.set_Y_name(starty);
	hashed_genotypes=hash_genotypes(file_buffer, startx, starty, l2o, false);
	gsl_vector *v=gsl_vector_alloc(7);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<double> unif(-0.01, 0.01);
    std::uniform_int_distribution<> unifint(0, 6);
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );

	//gsl_vector_set(v, unifint(eng), double(unif(eng))/20. );

    gsl_vector_set_all(v, 0 );
    double sum=rel_ll(v, &hashed_genotypes_vector);// << std::endl << std::endl;

    std::cerr << "Zeros:  " << sum << std::endl;

    if (around.size()==7){
        gsl_vector_set(v, 0, around[0] );
        gsl_vector_set(v, 1, around[1] );
        gsl_vector_set(v, 2, around[2] );
        gsl_vector_set(v, 3, around[3] );
        gsl_vector_set(v, 4, around[4] );
        gsl_vector_set(v, 5, around[5] );
        gsl_vector_set(v, 6, around[6] );
        std::cerr << around[0] << ", " << around[1] << ", " << around[2] << ", " << around[3] << ", " << around[4] << ", " << around[5] << ", " << around[6] << std::endl;
    } else {
        std::cerr << "NO! Size seven you git!\n";
    }
    sum=rel_ll(v, &hashed_genotypes_vector);// << std::endl << std::endl;
    std::cerr << "Test:  " << sum << std::endl;

    /*
    for (int x=0; x<101; x++){
	    gsl_vector_set(v, axis, around[axis]+float(x-50.)/500. );
        //std::cout << std::endl << std::endl;
        //for (int y=0; y<7; y++) std::cout << gsl_vector_get(v,y) << ", ";
        sum=rel_ll(v, &hashed_genotypes_vector);// << std::endl << std::endl;
        std::cerr << around[axis]+float(x-50.)/500. << ", " << sum << std::endl;
    }
    */
    
    for (int x=0; x<1000000; x++){
        axis=unifint(eng);
        for (int y=0; y<7; y++) 
        {
        gsl_vector_set(v, y, around[y]+unif(eng) );
        }
        double ll=rel_ll(v, &hashed_genotypes_vector);
        if (!std::isnan(ll) ){
            for (int y=0; y<7; y++) std::cout << gsl_vector_get(v,y) << ", ";
            std::cout << ll << std::endl;
        }
    }
         
    return 0;					//Since everything worked, return 0!.
}

#else

int surface(int argc, char *argv[]) {
    printf("Having linking problems with MPI...");
    return 0;
}
#endif
