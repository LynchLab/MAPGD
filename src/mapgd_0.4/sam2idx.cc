/* 

command filter:

*/

#include "sam2idx.h"

int sam2idx(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	bool binary=false;

	std::string headerfile="";

	env_t env;
	env.setname("mapgd genotype");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("convert a sam header file into an idx file.");
	env.optional_arg('H',"header",	&headerfile,	&arg_setstr, 	"an error occured", "sets the index file (required to use mpileup)");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'h', "help", &env, 		&flag_help, 	"an error occured while displaying the version message.", "prints the program version");

	if (parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	Flat_file <File_index> index_file;
        File_index index;

        if (headerfile.size()!=0) {
                std::fstream header;
                header.open(headerfile.c_str(),  std::fstream::in);
                index.from_sam_header( header );
        } else {
                index.from_sam_header( std::cin );
	}

	index_file.open(std::fstream::out);
	index_file.write_header(index);
	index_file.write(index);
	index_file.close();
	return 0;
}
