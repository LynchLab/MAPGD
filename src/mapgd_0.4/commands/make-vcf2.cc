/* 
*/

#include "make-vcf2.h"

#define BUFFER_SIZE 500



int make_vcf2(int argc, char *argv[])
{

#ifndef NOHTS
#ifndef NOLZ4
	/* All the variables that can be set from the command line */

	std::string statefile="";
	std::string namefile="";
	std::string indexname="";
	std::string outfile="";

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
	env.set_name("mapgd writevcf");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Convert state output into a vcf file (should be merged with make_vcf.");
	env.optional_arg('o',"output", 	outfile,"an error occurred while setting the name of the output file.", "the output vcf file (default stdout).");
	env.required_arg('s',"statefile", statefile,	"an error occurred while setting the name of the input file.", "the input 'stf' file (default none).");
	env.required_arg('n',"namefile", namefile,	"an error occurred while setting the name of the input file.", "the input 'stf' file (default none).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Flat_file <State> state_in;
	Flat_file <Sample_name> name_in;
	std::cerr << __LINE__ << std::endl;

	Vcf_file vcf_out;
	Sample_name names;

	state_in.open(statefile.c_str(), READ);
	name_in.open(namefile.c_str(), READ);

	State state;
	state=state_in.read_header();
	names=name_in.read_header();
	state_in.read(state);
	state_in.close();

	std::cerr << state.sample_size() << ", " << state.genome_size() << std::endl;

	name_in.read(names);

	size_t N = state.sample_size();
	size_t sites = state.genome_size();
	uint32_t *P1=new uint32_t [N];
	uint32_t *P2=new uint32_t [N];

	std::cerr << __LINE__ << std::endl;

	Vcf_data vcf;

	File_index index;
	index.add_id("1", 32*sites+1);

	vcf.set_header(index, names.sample_names);

	vcf_out.open(outfile.c_str(), WRITE);
	vcf_out.write_header(vcf);

	Allele allele;
	Population pop(names);

	std::cerr << __LINE__ << std::endl;
	static uint32_t mask[32]={0x00000001,   0x00000002,     0x00000004,     0x00000008,
		0x00000010,     0x00000020,     0x00000040,     0x00000080,
		0x00000100,     0x00000200,     0x00000400,     0x00000800,
		0x00001000,     0x00002000,     0x00004000,     0x00008000,
		0x00010000,     0x00020000,     0x00040000,     0x00080000,
		0x00100000,     0x00200000,     0x00400000,     0x00800000,
		0x01000000,     0x02000000,     0x04000000,     0x08000000,
		0x10000000,     0x20000000,     0x40000000,     0x80000000};

	state.rewind();

	allele.ref=0;
	allele.major=1;
	allele.minor=2;

	for (size_t x=0; x<sites; ++x)
	{
		state.uncompress(P1, P2);
		for (size_t k=0; k<32; ++k)
		{
			uint32_t count=0;
			for (size_t y=0; y < N; ++y) 
			{
				char g=( (P1[y] & mask[k]) !=0)+( (P2[y] & mask[k]) !=0);
				switch (g) 
				{
					case 0:
						pop.likelihoods[y].mm=1;
						pop.likelihoods[y].Mm=0;
						pop.likelihoods[y].MM=0;
					break;
					case 1:
						pop.likelihoods[y].mm=0;
						pop.likelihoods[y].Mm=1;
						pop.likelihoods[y].MM=0;
						count+=1;
					break;
					case 2:
						pop.likelihoods[y].mm=0;
						pop.likelihoods[y].Mm=0;
						pop.likelihoods[y].MM=1;
						count+=2;
					break;
				}
				pop.likelihoods[y].N=100;
			}
			allele.freq=1.-float(count)/float(2*N);
			allele.set_abs_pos(x*32+k+1);
			pop.set_abs_pos(32*x+k+1);
			
			vcf.put(index, allele, pop);
			vcf_out.write(vcf);
		}
	}

	vcf_out.close();
	return 0;					//Since everything worked, return 0!.
#else 
	std::cerr << "You must have lz4 to use this command. If you are using linux you can obtain lz4 by typing apt-get install liblz4-dev\n";	
	return 0;
#endif
#else 
	std::cerr << "You must have htslib to use this command. If you are using linux you can obtain htslib by typing apt-get install htslib-dev\n";	
	return 0;
#endif
}
