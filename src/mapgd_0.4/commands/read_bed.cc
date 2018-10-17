/* 
*/

#include "read_bed.h"

#define BUFFER_SIZE 500

int read_bed(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string bedfile="";
	std::string sttfile="";

	bool binary=false;
	int n=0;

	Environment env;
	env.set_name("mapgd readbed");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert a bed file into a stt file.");
	env.required_arg('n',"number",  n,	"an error occurred setting the number of individuals encoded.", "the number of individuals encoded in the bed file (default 0)");
	env.optional_arg('o',"output",  sttfile,	"an error occurred while setting the name of the output file.", "the out file for the program (default stdout)");
	//env.positional_arg('i',"input", bedfile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdin)");
	env.optional_arg('i',"input", bedfile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'b', "binary", 	&binary,	&flag_set, 	"an error occurred while setting binary flag.", "output in binary format");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Bed_file in_bed;
	Bed_data bed(n);

	in_bed.open(bedfile.c_str(), READ);
	
#ifndef NOLZ4
	State states(n);
	Flat_file <State> state_out;

	if (sttfile!="") state_out.open(sttfile.c_str(), binary ? WRITE | BINARY : WRITE);
	else state_out.open(binary ? WRITE | BINARY : WRITE);
	
	std::cerr << in_bed.table_is_open() << std::endl;

	size_t max_size=system_memory() >> 1;

	bool no_header=false;

	size_t junk=0;
	while(in_bed.table_is_open() )
	{
		in_bed.read(bed);
		if (in_bed.table_is_open() )  
		{
			bed.get(states);
		}
		if ((junk++)%10000==0)	std::cerr << states.buffer_size() << ", " << max_size << std::endl;
		if (states.buffer_size() > max_size) {
			if (!no_header)
			{
				no_header=true;
				states.set_streaming(true);
				state_out.write_header(states);
			}
			states.cache();
			state_out.write(states);
		}
	}
	if (!states.cached() )
	{
		//states.finalize();
		states.cache();
		if (!no_header) 
		{
			state_out.write_header(states);
		}
		std::cerr << "not cached? \n";
		state_out.write(states);
	} else {
		std::cerr << "cached!?!? \n";
	}

	/*states.finalize();
	states.cache();*/
	state_out.close();

	in_bed.close();
#endif
	return 0;					//Since everything worked, return 0!.
}
