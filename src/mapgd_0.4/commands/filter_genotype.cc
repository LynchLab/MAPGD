/* 

command filter:

*/

#include "filter_genotype.h"

int filter_genotype(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string in_file="", out_file="";
	std::vector <size_t>  mins, maxs;
	std::vector <size_t>::iterator  min_it, max_it;

	int E;

	bool binary=false, twofold=false;

    int c, C;
    double min_maf=0, max_maf=1;

	Environment env;
	env.set_name("mapgd filter");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Filter sites in '.gcf' files based on criteria.");
	
	env.optional_arg('i',"input", 	in_file, "please provide a filename.", "the name of the input file (default stdin).");
	env.optional_arg('o',"output", 	out_file, "please provide a filename.", "the name of an output file (default stdout).");
	env.optional_arg('q',"minmaf", 	min_maf, 	"please provide an float.", "minimum minor allele frequency (default 0.0).");
	env.optional_arg('Q',"minmaf", 	max_maf, 	"please provide an float.", "maximum minor allele frequency (default 0.5).");
	env.optional_arg(0,"mincovs", 	mins, 	"please provide a list of integers.", "a list of minimum individual coverage (default none).");
	env.optional_arg(0,"maxcovs", 	maxs, 	"please provide a list of integers.", "a list of maximum individual coverage (default none).");
	env.optional_arg('c',"mincov", 	c, 	"please provide an integer.", "minimum individual coverage (default none).");
	env.optional_arg('C',"maxcov", 	C, 	"please provide an integer.", "maximum individual coverage (default none).");
	env.optional_arg('E',"exact", 	E, 	"please provide an integer.", "coverage of exactly N.");

	env.flag(	'2', "twofold",	&twofold,	&flag_set, 	"none.", "two fold coverage cutoff.");
	env.flag(	'b', "binary", 	&binary,	&flag_set, 	"none.", "output in binary mode (fast).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Population s;
	Genotype g, clear;
	Indexed_file <Population> gcf_in, gcf_out;

	if (in_file.size()==0)	gcf_in.open(std::fstream::in);
	else gcf_in.open(in_file.c_str(), std::fstream::in);

	if (out_file.size()==0)	gcf_out.open( binary ? WRITE |  BINARY : WRITE );
	else gcf_out.open(out_file.c_str(), binary ? WRITE | BINARY : WRITE );


	s=gcf_in.read_header();

	gcf_out.set_index(gcf_in.get_index() );
	gcf_out.write_header(s);

	gcf_in.read(s);

        std::vector <Genotype>::iterator it, begin, end;             //!< Genotypic likelihood

	while( gcf_in.table_is_open() ){
		end=s.likelihoods.end();
		min_it=mins.begin();
		max_it=maxs.begin();
		for (it=s.likelihoods.begin(); it!=end; ++it)
		{
            //if ?
    		//	if(it->N<*c) *it=clear;
	    	//	else if(it->N>*C) *it=clear;
    	    if(it->N<c) *it=clear;
	    	else if(it->N>C) *it=clear;

			//++c;
			//++C;
	    }	
        if (s.m <= max_maf && s.m >= min_maf )
    		gcf_out.write(s);
		gcf_in.read(s);
	}
	gcf_out.close();
	return 0;					//Since everything worked, return 0!.
}
