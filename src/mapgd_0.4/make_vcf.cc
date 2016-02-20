/* 
*/

#include "vcf.h"

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
	env.set_name("mapgd allele");
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

	Vcf_file vcf_out;
	
	Allele allele_in;
	Locus locus_in;

	Vcf_header 
	Vcf_
	if (infile.size()!=0) {					//Iff a filename has been set for infile
		pro_in.open(infile.c_str(), std::fstream::in);	
		if (!pro_in.is_open() ) {				//try to open a profile of that name.
			print_usage(env);			//Print help message on failure.
		} 
		locus_in=pro_in.read_header();
	}
	else {
		pro_in.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
		locus_in=pro_in.read_header();
	};
	if (!pro_in.table_is_open() ) 
	{
		fprintf(stderr, "%s:%d. Error: cannot open file.\n",__FILE__, __LINE__);
		exit(0);
	}


	bool binary=false;

	if (outfilepro.size()!=0) {				//Same sort of stuff for the outf
		if (binary) {
			pro_out.open(outfilepro.c_str(), std::fstream::out | std::fstream::binary);
			if (!pro_out.is_open() ) {print_usage(env); exit(0);}
		} else {
			pro_out.open(outfilepro.c_str(), std::fstream::out);
			if (!pro_out.is_open() ) {print_usage(env); exit(0);}
		}
		locus_out.set_sample_names(locus_in.get_sample_names() );
		pro_out.set_index(pro_in.get_index() );
		pro_out.write_header(locus_out);
	};


	/* this is the basic header of our outfile, should probably be moved over to a method in Allele.*/
	locus_in.maskall();						//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<locus_in.get_sample_names().size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	std::vector <float_t> sum_gofs(ind.size() );
	std::vector <float_t> gofs_read(ind.size() );

	models model;

	if (outfile.size()!=0) {
		map_out.open(outfile.c_str(), std::ios::out);
		if (!map_out.is_open() ) print_usage(env);
	} else 	map_out.open(std::fstream::out);

	Allele buffer_mle[BUFFER_SIZE]; 
	Locus buffer_locus[BUFFER_SIZE];
	std::fill_n(buffer_locus, BUFFER_SIZE, locus_in);

	map_out.set_index(pro_in.get_index() );
	map_out.write_header(buffer_mle[0]);

	uint32_t all_read=0;

	while (true){			//reads the next line of the pro file. pro.read() retuerns 0
		uint32_t c=0, readed=0;
		bool estimate_me=true;
		#ifndef NOOMP
		#pragma omp parallel private(c, model, estimate_me) 
		#endif
		{
			#ifndef NOOMP
			#pragma omp for
			#endif
			for (uint32_t x=0; x<BUFFER_SIZE; ++x){
				#ifndef NOOMP
				#pragma omp critical
				#endif
				{
					c=readed;		//Turn on the ability to read data from all clones in 
					if(pro_in.read(buffer_locus[c]).table_is_open() ){
						readed++;	//reads the next line of the pro file. pro.read() retuerns 0
						estimate_me=1;
					}
					else estimate_me=0;
				}
				if(estimate_me) {
					std::vector <float_t> gofs(ind.size() );
		//			buffer_locus[c].unmaskall();
					buffer_locus[c].unmask(ind);

					buffer_mle[c]=estimate (buffer_locus[c], model, gofs, MIN, EMLMIN, MINGOF, MAXPITCH, newton);
					#ifndef NOOMP
					#pragma omp critical
					#endif
					if (2*(buffer_mle[c].ll-buffer_mle[c].monoll)>=22){
						for (size_t i=0; i<sum_gofs.size(); i++){
							sum_gofs[i]+=gofs[i];
							if (gofs[i]!=0) gofs_read[i]++;
						}
					}
				}
			}

		}
                for (uint32_t x=0; x<readed; ++x){
			map_out.write(buffer_mle[x]);
			if (buffer_mle[x].gof<-MINGOF) buffer_locus[x].maskall(); 
			if (pro_out.is_open() ){
				pro_out.write(buffer_locus[x]);
			}
		}
		if (readed!=BUFFER_SIZE){break;}
		all_read+=readed;
		if (all_read>stop){break;}
	}
	map_out.close_table();
	if ( not(noheader) ) {
		Flat_file <Sample_gof> gof_file;
		gof_file.open_from(map_out);
		gof_file.write_header(Sample_gof() );
		for (size_t x=0; x<ind.size(); ++x) gof_file.write(Sample_gof(locus_in.get_sample_names()[ind[x]], sum_gofs[x]/(float_t(gofs_read[x]) ) ) );
		gof_file.close();
	}
	pro_in.close();
	if (pro_out.is_open()) pro_out.close();		//Closes pro_out iff pro_out is open.
	return 0;					//Since everything worked, return 0!.
}
