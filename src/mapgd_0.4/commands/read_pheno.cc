#include "read_pheno.h"

int read_pheno(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string in_file="";
	std::string out_file="";

	bool binary=false;
	bool state=false;

	Environment env;
	env.set_name("mapgd readpheno");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert a pheno file into a phe file.");
	env.optional_arg('o',"output", 	out_file,"an error occurred while setting the name of the output file.", "the output file for the program (default stdout)");
	env.positional_arg('i',"input", in_file,	"an error occurred while setting the name of the input file.", "the input file for the program (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Flat_file <Phenotype> pheno_out;
	
	if (out_file!="") pheno_out.open(out_file.c_str(),  binary ? WRITE | BINARY : WRITE);
	else pheno_out.open(WRITE);

	Plink_file plink_file;
	if (in_file!="") plink_file.open(in_file.c_str(), READ);
	else plink_file.open(READ);

	Plink_data plink;
	plink=plink_file.read_header();
	
	Phenotype pheno(plink.trait_size() );
	

	while(true)
	{
		plink_file.read(plink);
		if (plink_file.table_is_open() ) plink.get(pheno);
		else break;
	}

	pheno_out.write_header(pheno);
	pheno_out.write(pheno);
	plink_file.close();
	pheno_out.close();
	return 0;
}
