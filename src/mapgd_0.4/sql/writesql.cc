/* 

command filter:

*/

#include "writesql.h"

int writesql(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	env_t env;
	env.setname("mapgd writesql");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("writes mapgd output to an sql database.");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'h', "help", &env, 		&flag_help, 	"an error occured while displaying the version message.", "prints the program version");

	if (parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

		
	return 0;					//Since everything worked, return 0!.
}
