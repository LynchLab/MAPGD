/* 

command filter:

*/

#include "map2genotype.h"

//This should probably be changed over to a likelihood model class
genotype baysian_genotype(const int &major, const int &minor, const float_t &freq, const float_t &error, const quartet_t &quart)
{
	float_t lMM, lMm, lmm, N;

	N=count(quart);

	count_t M=quart.base[major], m=quart.base[minor], E=N-M-m;

	float_t ln_homozygous_correct=log(1.-error);
	float_t ln_heterozygous_correct=log( (1.-error)/2.+error/6.);
	float_t not_correct=log(error/3.);

	lMM=M*ln_homozygous_correct+not_correct*(m+E);
	lmm=(M+m)*ln_heterozygous_correct+not_correct*(E);
	lmm=m*ln_homozygous_correct+not_correct*(M+E);

	float_t norm=log(exp(lMM)+exp(lMm)+exp(lmm) );

	lMM=-lMM+norm;
	lMm=-lMm+norm;
	lmm=-lmm+norm;

	return genotype(lMM, lMm, lmm, N);	
}

void get_genotypes(const Allele &allele, const Locus &locus, population_genotypes &genotypes )
{
	std::vector <genotype>::iterator l_it=genotypes.likelihoods.begin();
	std::vector <quartet_t>::const_iterator q_it=locus.cbegin(), q_end=locus.cend();
	while (q_it != q_end){
		*l_it=baysian_genotype(allele.major, allele.minor, allele.freq, allele.error, *q_it);
		q_it++;
		l_it++;
	}
	genotypes.m=1.-allele.freq;
	genotypes.major=allele.major;
	genotypes.minor=allele.minor;
	genotypes.id0=allele.id0;
	genotypes.id1=allele.id1;
} 
 

int map2genotype(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	bool binary=false;

	std::string mapname="", proname="";

	env_t env;
	env.setname("mapgd genotype");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("convert a map and pro file to an individual genotype file.");
	env.required_arg('m',"map", &mapname, 	&arg_setstr, "please provide a float.", "the name of the map file.");
	env.required_arg('p',"pro", &proname, 	&arg_setstr, "please provide a float.", "tne name of the pro file.");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'h', "help", &env, 		&flag_help, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'b', "binary", &binary, 	&flag_set, 	"an error occured while displaying the version message.", "binary output");

	if (parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	/* The files we will need*/
	Indexed_file <Allele> map_in;
	Indexed_file <Locus> pro_in;
	Indexed_file <population_genotypes> gcf_out;

	/* Open input files based on file name*/
	map_in.open(mapname.c_str(), std::ios::in);
	pro_in.open(proname.c_str(), std::ios::in);

	/* Just open the output file to std::cout*/
	gcf_out.open(binary ? std::ios::out | std::ios::binary : std::ios::out);

	/* We need one instance of each of the classes for reading and 
	 * writing to files */
	Locus pro_record;
	Allele map_record;
	population_genotypes gcf_record;

	/* Read the headers of the files */
	map_record=map_in.read_header();	
	pro_record=pro_in.read_header();	

	/* Set the sample names for the gcf file from the sample names in the pro_file*/
	gcf_record.set_sample_names(pro_record.get_sample_names() );
	gcf_out.set_index(map_in.get_index() );

	/* Write the header */
	gcf_out.write_header(gcf_record);

	id1_t map_pos;	
	map_in.read(map_record);
	
	while(map_in.table_is_open() ){
	/* a read/write cycle */
		map_pos=map_in.get_pos(map_record);
		while(pro_in.get_pos(pro_record)<map_pos && !pro_in.eof() ){
			pro_in.read(pro_record);
//			std::cerr << pro_record << std::endl;
		}
		if (map_pos==pro_in.get_pos(pro_record) ){
			get_genotypes(map_record, pro_record, gcf_record);
			gcf_out.write(gcf_record);
		}
		map_in.read(map_record);
	}
	gcf_out.close();

	return 0;					//Since everything worked, return 0!.
}
