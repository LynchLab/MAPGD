// Updated on 01/19/16

#include "linkage_disequilibrium.h"

/*

Program linkage_disequilibrium.cc to estimate linkage disequilibrium (LD) between 
polymorphic sites from high-throughput sequencing data of multiple 
diploid individuals from a population by a maximum-likelihood (ML) method.
Random union of gametes is assumed in the estimation.  In addition,
all sequence reads are assumed to independently cover at most one
of the two polymorphic sites of interest.

*/

//TODO some performance profiling on buffer sizes.
#define BUFFER_SIZE 1000

using namespace std;

count_t
count_sites(const Population &X, const Population &Y)
{
	count_t ret=0;
	std::vector<Genotype>::const_iterator it_X=X.liklihoods.cbegin();
	std::vector<Genotype>::const_iterator it_Y=Y.likelihoods.cbegin();
	std::vector<Genotype>::const_iterator end=Y.likelihoods.cend();
	while (it_Y!=end){
		ret+=(count(*it_X).N!=0) & (count(*it_Y).N!=0);
		++it_X;
		++it_Y;
	}
	return ret;
}

float_t guesstimate_D (const Population &P1, const Population &P2)
{
	return 0;	
}

Linkage estimate_D (const Population &P1, const Population &P2, float_t &D)
{
	float_t R=H0(P1, P2, D);
	while (R>0.0001){
		R=H0(P1, P2, D);
		D-=R/J00(P1, P2, D)
		
        }
	
	Linkage est;
        est.set_abs_pos(P1.get_abs_pos() );
        est.set_abs_pos_y(P2.get_abs_pos() );
        est.set_D(D);
        est.set_fit(lnL(P1, P2, D) );
        est.set_null(lnL(P1, P2, 0 ) );

        return(est);
}
	
int linkage_disequilibrium(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string pro_name="";
	std::string map_name="";

	std::vector <int> ind;

	int max_d = INT_MAX;
	int min_dist=0;
	double min_number = 10.0;
	int print_help = 0;

	Environment env;
	env.set_name("mapgd linkage");
	env.set_version(VERSION);
	env.set_author("Takahiro Maruki and Matthew Ackerman");
	env.set_description("Uses a maximum likelihood approach to estimate gametic phase disequilibrium from population data.");

	env.optional_arg('i',"in", 	gcf_name,	"please enter a string.", "the input 'gcf' file (default stdin).");
	env.optional_arg('o',"out", 	gcf_name,	"please enter a string.", "the input 'gcf' file (default stdin).");
	env.optional_arg('M',"min_n", 	min_number, 	"please enter a number.", "the minimum number of individuals at a site need to calculate LD (default: 10).");
	env.optional_arg('D',"max_d", 	max_d,		"please enter a number.", "the maximum distance between sites for LD calculation (default: none).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.
	
	// Default values of the options

	Indexed_file <Population> gcf_in;

	if gcf_name
	gcf_in.open(map_name.c_str() ,std::ios::in);	// Try to open the input file
	else 	

	Indexed_file <Linkage> ld_out;	
	ld_out.open(std::ios::out);

	Population allele1=gcf_in.read_header();
	Population allele2=allele1;	
	File_index index=gcf_in.get_index();
	ld_out.set_index(index);
	Linkage linkage;
	ld_out.write_header(linkage);

	std::list <Population> allele_list;

	Population allele_buffer1[BUFFER_SIZE];
	Population allele_buffer2[BUFFER_SIZE];

	std::fill_n(allele_buffer1, BUFFER_SIZE, allele1);
	std::fill_n(allele_buffer2, BUFFER_SIZE, allele1);

	Linkage linkage_buffer[BUFFER_SIZE]; 

	std::fill_n(linkage_buffer, BUFFER_SIZE, linkage);

	count_t nsample=allele1.get_sample_names().size();

	id1_t loc=0, read=0;

	allele_list.insert(allele_list.end(), BUFFER_SIZE, allele1);

	std::list <Population>::iterator s_allele=allele_list.begin(), e_allele=allele_list.begin(), end_allele=allele_list.end();

	while(e_allele!=end_allele){
		gcf_in.read(*(e_allele) );
		e_allele++;
	}

	e_allele=allele_list.begin();
						
	do {
		read=0;
		do {
			size_t number;
			id1_t pos1=allele1.get_abs_pos();
			id1_t pos2=allele2.get_abs_pos();
			id0_t scf1=index.get_id0(pos1);
			id0_t scf2=index.get_id0(pos2);

			//TODO Fix my terrible code! This code block checks to see if the two sites (locus1 and locus2) pass critera for calculating LD.
			if ( (pos2-pos1)<max_d && scf1==scf2) {
				if ( (pos2-pos1)>min_dist ) {
					number=count_sites(allele1, allele2);
					//std::cerr << number << ", " << pos1 << ", " << pos2 << ", " << pos2-pos1 << " ding " << std::distance(e_allele,end_allele) << " -- " << std::distance(s_allele, end_allele) <<  ".\n";
					if ( number >= min_number ) {
						allele_buffer1[read]=allele1;
						allele_buffer2[read]=allele2;
						read++;
						//size_t Ni=count_sites(locus_buffer1[x], locus_buffer2[x]);
					}
				}
			} else { 
				//If the distance between the two sites is too great, move the first site up.
				s_allele++;
				if (s_allele!=end_allele && e_allele!=end_allele){
					allele1=*s_allele;
					e_allele=s_allele;
					allele_list.pop_front();
				} else {
					e_allele=end_allele;
				}
			}

			//Make sure site two isn't at the end of our buffer.
			if (e_allele!=end_allele) {
				//If it isn't, just move site two up..
				allele2=*(e_allele);
				e_allele++;
			} else if (gcf_in.table_is_open() ) {
				//If it is we have to increase the size of our buffer.
				//std::cerr << "Buffer read\n";
				std::list <Population>::iterator new_allele=e_allele;
				new_allele--;
				allele_list.insert(e_allele, BUFFER_SIZE, allele1);

				end_allele=allele_list.end();

				e_allele=new_allele;

				while(new_allele!=end_allele){
					gcf_in.read( *(new_allele) );
					id1_t map_pos=gcf_in.get_pos(*new_allele);
					new_allele++;
				}

				++e_allele;

				allele2=*(e_allele);
			} else {
				//std::cerr << "Moving locus1\n";
				s_allele++;
				if (s_allele!=end_allele){
					allele1=*s_allele;
					e_allele=s_allele;
					allele_list.pop_front();
				} else {
					e_allele=end_allele;
				}
			}
		} while (read<BUFFER_SIZE && (gcf_in.table_is_open() || s_allele!=end_allele) );
	
		#pragma omp for
		for (uint32_t x=0; x<BUFFER_SIZE; ++x){
			if (x<read){
				//size_t Ni=count_sites(allele_buffer1[x], allele_buffer2[x]);
				//if (Ni>
				linkage_buffer[x] = estimate_D( allele_buffer1[x], allele_buffer2[x] ,guesstimate_D(allele_buffer1[x], allele_buffer2[x]) );
			}
		}
		for (size_t c=0; c<read; ++c){ 
			ld_out.write(linkage_buffer[c]);
		}
	} while (gcf_in.table_is_open() || s_allele!=end_allele);
	ld_out.close();
	return 0;
}
