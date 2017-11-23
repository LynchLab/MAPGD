/* 
*/

#include "read_vcf.h"

#define BUFFER_SIZE 500

int read_vcf(int argc, char *argv[])
{

#ifndef NOHTS
	/* All the variables that can be set from the command line */

	std::string vcffile="";
	std::string gcffile="";
	std::string indexfile="";

	bool binary=false;
	bool state=false;

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
	env.set_name("mapgd readvcf");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert a vcf file into a gcf or stt file.");
	env.optional_arg('o',"output", 	gcffile,"an error occurred while setting the name of the output file.", "the output file for the program (default stdout)");
	env.optional_arg('H',"header", 	indexfile,"an error occurred while setting the name of the output file.", "sets the index file");
	env.positional_arg('i',"input", vcffile,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	's', "state", 	&state,		&flag_set, 	"an error occurred while displaying the help message.", "output state file");
	env.flag(	'b', "binary", 	&binary,	&flag_set, 	"an error occurred while displaying the help message.", "output in binary format");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.


	Vcf_file vcf_in;
	Vcf_data vcf;
	File_index index;

	vcf_in.open(vcffile.c_str(), READ);
	vcf=vcf_in.read_header();
	
	if (indexfile=="") {
		index=vcf.get_index();
		std::cerr << "This vcf contains " << index.get_names_size() << " scaffold names. If this is not correct try specifying an index file (-H) created with the sam2idx command." << std::endl;
		if (index.get_names_size()==0) exit(0);
	}

	std::vector<std::string> names=vcf.get_sample_names();


	if (!state){	
		Indexed_file <Population> gcf_out;
		Indexed_file <Allele> map_out;
	
		std::vector <std::string> header_line={"@ID0","ID1","REF"};
		header_line.insert(header_line.end(), names.begin(), names.end() );	

		Population pop(header_line);
		if (gcffile!="") gcf_out.open(gcffile.c_str(),  binary ? WRITE | BINARY : WRITE);
		else gcf_out.open(WRITE);
		gcf_out.set_index(index);
		gcf_out.write_header(pop);

		Allele allele(header_line);

//		map_out.open(mapfile.c_str(), WRITE);
//		map_out.write_header(allele);

		while(vcf_in.table_is_open() )
		{
			vcf_in.read(vcf);
//			vcf.get(index, pop);
//			gcf_out.write(pop);
		}

		vcf_in.close();
		gcf_out.close();
	} else {
#ifndef NOLZ4
		State states(vcf.get_sample_names().size() );
		Flat_file <State> state_out;

		if (gcffile!="") state_out.open(gcffile.c_str(), binary ? WRITE | BINARY : WRITE);
		else state_out.open(binary ? WRITE | BINARY : WRITE);

		while(vcf_in.table_is_open() )
		{
			vcf_in.read(vcf);
			if (vcf_in.table_is_open() ) vcf.get(states);
	//		state_out.write(states);
		}
		states.finalize();
		states.cache();
		state_out.write_header(states);
		state_out.write(states);
		state_out.close();

		vcf_in.close();
#else 
		std::cerr << "Whoops! No LZ4 compression...\n";
#endif
	}
	return 0;					//Since everything worked, return 0!.
#else 
	std::cerr << "You must have htslib to use this command. If you are using linux you can obtain htslib by typing apt-get install htslib-dev\n";	
	return 0;
#endif
}
