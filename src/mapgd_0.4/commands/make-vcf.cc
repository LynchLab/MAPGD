/* 
*/

#include "make-vcf.h"

#define BUFFER_SIZE 500

int vcf(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string infile="";
	std::string outfile="";
	std::string outfilepro;
	std::string indexname="";

	bool verbose=false;
	bool ldformat=false;
	bool quite=false;
	bool noheader=false;
	bool newton=false;

	int rnseed=3;

	float_t EMLMIN=0.001;
	count_t MIN=4;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	id1_t skip=0;
	id1_t stop=UINT32_MAX;

	std::vector <size_t> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */
	//std::cerr "If you are using .... please cite ...."

	Environment env;
	env.set_name("mapgd vcf");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert mapgd output into a vcf file.");
	env.optional_arg('o',"output", 	outfile,"an error occurred while setting the name of the output file.", "the output file for the program (default stdin).");
	env.positional_arg('i',"input", infile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdout).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Allele> map_in;
	Indexed_file <Locus> pro_in;

//	Vcf_file vcf_out;
	
	Allele allele_in;
	Locus locus_in;

	return 0;					//Since everything worked, return 0!.
}
