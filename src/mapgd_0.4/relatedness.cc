/* 

Program estimateIndcpp:
*/

#include "relatedness.h"

#define BUFFER_SIZE 500
#define PRAGMA

int estimateRel(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string infile="";
	std::string outfile="";
	bool verbose=false;
	bool ldformat=false;
	bool quite=false;
	bool noheader=false;
	float_t EMLMIN=0.001;
	count_t MIN=4;
	float_t A=0.00;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	count_t skip=0;
	count_t stop=-1;

	std::vector <size_t> ind;

	env_t env;
	env.setname("mapgd relatedness");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	profile pro, pro_out;		//profile is a fairly complete class that lets us read and write from pro files, 
					//which are files containing set of read 'quartets' that specify the number of 
					//A,C,G and T read at some specific location in a genome. See proFile.h for more info.

	//gcfile out;
	std::ostream *out=&std::cout;
	std::ofstream outFile;

	return 0;					//Since everything worked, return 0!.
}
