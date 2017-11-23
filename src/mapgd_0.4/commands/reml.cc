#include "reml.h"


int reml(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string relfile="";
	std::string traitfile="";
	std::string outfile="";
	std::string indexname="";

	bool verbose=false;
	bool quite=false;
	bool noheader=false;
	bool newton=false;

	int rnseed=3;

	Environment env;
	env.set_name("mapgd reml");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses a restricted maximum likelihood approach to estimate population level components of variation (does not currently work).");
	env.optional_arg('r',"seed", 	rnseed,	"please provide a number.", "random number seed (3).");

	env.optional_arg('o',"output", 	outfile,	"an error occurred while setting the name of the output file.", "the output file for the program (default stdin).");
	env.positional_arg('r',"relatedness",	relfile,	"No input file specified", "the input file for the program (default stdout).");
	env.positional_arg('t',"traits",	traitfile,	"No input file specified", "the input file for the program (default stdout).");

	env.flag(	'n',"newton", 	&newton,	&flag_set, 	"takes no argument", "use Newton-Raphson likelihood maximization (not working).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");
	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occurred while enabling verbose execution.", "prints more information while the command is running.");
	env.flag(	'q', "quiet", 	&quite,		&flag_set, 	"an error occurred while enabling quite execution.", "prints less information while the command is running.");


	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Flat_file <Relatedness> rel_in;
	//Flat_file <Trait> trait_in;

	Relatedness rel;
	//Trait trait;	

	if (relfile.size()!=0) {					//Iff a filename has been set for infile
		rel_in.open(relfile.c_str(), std::fstream::in);	
		if (!rel_in.is_open() ) {				//try to open a profile of that name.
			print_usage(env);			//Print help message on failure.
		} 
		rel=rel_in.read_header();
	}
	else {
		rel_in.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
		rel=rel_in.read_header();
	};
	
	/*if (traitfile.size()!=0) {					//Iff a filename has been set for infile
		triat_in.open(triatfile.c_str(), std::fstream::in);	
		if (!trait_in.is_open() ) {				//try to open a profile of that name.
			print_usage(env);			//Print help message on failure.
		} 
		trait=trait_in.read_header();
	}
	else {
		trait_in.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
		trait=triat_in.read_header();
	};*/

	bool binary=false;

	rel_in.close();
	return 0;					//Since everything worked, return 0!.
}
