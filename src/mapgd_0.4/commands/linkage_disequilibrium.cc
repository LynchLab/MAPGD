// Updated on 01/19/16

#include "linkage_disequilibrium.h"
#include <algorithm>
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
	std::vector<Genotype>::const_iterator it_X=X.likelihoods.cbegin();
	std::vector<Genotype>::const_iterator it_Y=Y.likelihoods.cbegin();
	std::vector<Genotype>::const_iterator end=Y.likelihoods.cend();
	while (it_Y!=end){
		ret+=(it_X->N!=0 and it_Y->N!=0);
		++it_X;
		++it_Y;
	}
	return ret;
}

float_t guesstimate_D (const Population &P1, const Population &P2)
{
	return 0.00;	
}

Linkage estimate_D (const Population &P1, const Population &P2, const float_t &cD)//, std::default_random_engine &re)
{
	float_t q1=P1.m;
	float_t q2=P2.m;
	float_t p1=1.-q1;
	float_t p2=1.-q2;
	double Dmin=-min(p1*p2, q1*q2);
	double Dmax=min(q1*p2, q2*p1);
	float_t D=cD;
	//float_t Dp=cD, Dpp=cD+(cD-Dmax)/5.;
	//float_t Hp=H0(P1, P2, Dp);
	//float_t Hpp=H0(P1, P2, Dpp);
	float_t R=H0(P1, P2, D);
	float_t J;
//	std::cerr << " Called with " << cD << std::endl;
	//D=(Dpp*Hp-Dp*Hpp)/(Hp-Hpp);
	//We're going to try the method of the Secant...
	int iter=0;
	while (fabs(R)>0.05 and iter < 20){
		R=H0(P1, P2, D);
		J=J00(P1, P2, D);
		/*if ( fabs(J) > 100000 ){
			if (J > 0 ) J=100000;
			else J=-100000;
		}*/
		D-=R/J;
//		std::cerr << D << " = " << R << ", " << J << ", " << J00(P1, P2, D) << std::endl;
		iter++;
/*		Hpp=Hp;
		Dpp=Dp;
		Dp=D;
		Hp=H0(P1, P2, Dp);
		D=(Dpp*Hp-Dp*Hpp)/(Hp-Hpp);
		std::cerr << D << " = " << Hp << ", " << (Dpp*Hp-Dp*Hpp)/(Hp-Hpp) << std::endl;
*/ 
	}
	if (D < Dmin or D > Dmax or iter == 20) {
		float min_lnL=lnL_NR(P1, P2, Dmin);
		float max_lnL=lnL_NR(P1, P2, Dmax);
		if ( not isnan(min_lnL) ){
			if (min_lnL < max_lnL) {
				D=Dmax;
			}
			else D=Dmin;
		} else D=Dmax;
	} 
//	std::cerr << P1.m << "/" << P2.m << ", " << Dmin << ", " << Dmax << ", " << D << std::endl;
	
/*	float_t scan_max=-FLT_MAX, scan_D=0;	
	for (int x=0; x<100; x++)
	{
		//std::cerr << Dmin << ", " << Dmax << ", " << Dmin+x/100.*(Dmax-Dmin) << ", " << lnL_NR(P1, P2, Dmin+x/100.*(Dmax-Dmin) ) << std::endl;
		if (lnL_NR(P1, P2, Dmin+x/100.*(Dmax-Dmin) ) > scan_max) 
		{
			scan_max=lnL_NR(P1, P2, Dmin+x/100.*(Dmax-Dmin) );
			scan_D=Dmin+x/100.*(Dmax-Dmin);
		}
	}
	std::cerr << "converged to " << D << "=" << lnL_NR(P1, P2, D) << " and scan max is " << scan_max << " at " << scan_D  << std::endl;
*/
	
	Linkage est;

        est.set_Ni( count_sites(P1, P2) );
        est.set_p(P1.m);
        est.set_q(P2.m);

        est.set_abs_pos(P1.get_abs_pos() );
        est.set_abs_pos_y(P2.get_abs_pos() );
        est.set_D(D);
        est.set_fit(lnL_NR(P1, P2, D) );
        est.set_null(lnL_NR(P1, P2, 0 ) );

        return(est);
}

#define TRACE_LINKAGE
	
int linkage_disequilibrium(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string gcf_name="";

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

	if (gcf_name.size()>0 ) gcf_in.open(gcf_name.c_str() ,std::ios::in);	// Try to open the input file
	else gcf_in.open(std::ios::in);	

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
	
//	std::default_random_engine re;	
					
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
					if ( number >= min_number ) {
						allele_buffer1[read]=allele1;
						allele_buffer2[read]=allele2;
						read++;
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
				//Moving locus1 ...
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
				linkage_buffer[x] = estimate_D( allele_buffer1[x], allele_buffer2[x] ,guesstimate_D(allele_buffer1[x], allele_buffer2[x]) );//, re );
			}
		}
		for (size_t c=0; c<read; ++c){ 
			ld_out.write(linkage_buffer[c]);
		}
	} while (gcf_in.table_is_open() || s_allele!=end_allele);
	ld_out.close();
	return 0;
}
