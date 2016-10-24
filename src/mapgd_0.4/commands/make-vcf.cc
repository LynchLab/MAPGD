/* 
*/

#include "make-vcf.h"

#define BUFFER_SIZE 500

int make_vcf(int argc, char *argv[])
{

#ifndef NOHTS
	/* All the variables that can be set from the command line */

	std::string gcffile="";
	std::string mapfile="";
	std::string outfile="";
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
	env.positional_arg('g',"gcffile", gcffile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdout).");
	env.positional_arg('m',"mapfile", mapfile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdout).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in;
	Indexed_file <Allele> map_in;

	Vcf_file vcf_out;

	gcf_in.open(gcffile.c_str(), READ);
	map_in.open(mapfile.c_str(), READ);

	Population pop=gcf_in.read_header();
	Allele allele=map_in.read_header();

	File_index index=gcf_in.get_index();

	Vcf_data vcf;

	vcf.set_header(gcf_in.get_index(), pop.get_sample_names() );

	vcf_out.open(outfile.c_str(), WRITE);
	vcf_out.write_header(vcf);

	while(map_in.table_is_open() )
	{
		gcf_in.read(pop);
		map_in.read(allele);
		vcf.put(index, allele, pop);
		vcf_out.write(vcf);
	}

	vcf_out.close();

	gcf_in.close();
	map_in.close();

	return 0;					//Since everything worked, return 0!.
#else 
	std::cerr << "You must have htslib to use this command. If you are using linux you can obtain htslib by typing apt-get install htslib-dev\n";	
	return 0;
#endif
}
