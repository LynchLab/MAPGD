#include "convert.h"
#include <iostream>


//BY DEFAULT THIS SHOULD DO THE IO STREAM!!!
int convert(int argc, char *argv[])
{

        /* All the variables that can be set from the command line */

        std::string filename1="";
        std::string filename2="";

        env_t env;
        env.set_name("mapgd convert");
        env.set_version(VERSION);
        env.set_author("Matthew Ackerman");
        env.set_author("A small utility to convert a text file to a binary file.");
	env.required_arg('i',"input", 	&filename1,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdout).");
	env.required_arg('o',"input", 	&filename2,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdout).");
	
	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	streamtable(filename1.c_str(), filename2.c_str());
	return 0;
};
