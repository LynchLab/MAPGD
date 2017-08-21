/* 
*/

#include "simulate.h"

int simulate(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string vcffile="";
	std::string gcffile="";
	std::string indexfile="";

	Environment env;
	env.set_name("mapgd simulate");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert mapgd output into a vcf file.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Flat_file <State> state_file;
	state_file.read_header();
//	state_file.open(WRITE | BINARY );
	State state1(10), state2(10);

	uint32_t state1A[10]={10,9,8,7,6,5,4,3,2,1};
	uint32_t state1B[10]={1,2,3,4,5,6,7,8,9,10};
	uint32_t state2A[10]={0};
	uint32_t state2B[10]={0};

	state1.compress(state1A, state1B);
	state1.compress(state1B, state1A);

	state1.cache();

/*	state2.compress(state1A, state1B);
	state2.compress(state1B, state1A);
*/

	//swap(state1, state2);
	state2.cache();

	state2.uncompress(state2A, state2B);
	std::cerr << state2A[0] << ", " << state2B[9] << std::endl;
	state2.uncompress(state2A, state2B);
	std::cerr << state2A[0] << ", " << state2B[9] << std::endl;

	state2.rewind();

	state2.uncompress(state2A, state2B);
	std::cerr << state2A[0] << ", " << state2B[9] << std::endl;
	state2.uncompress(state2A, state2B);
	std::cerr << state2A[0] << ", " << state2B[9] << std::endl;

	state2.rewind();

	state2=state1;	
//	state_file.write_header(state2);
//	state_file.write(state2);
//	state_file.close();

	return 0;
}
