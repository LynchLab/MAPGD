// Updated on 01/19/16

#include "PopLD.h"
#include <iterator>     // std::distance

/*

Program PopLD.cc to estimate linkage disequilibrium (LD) between 
polymorphic sites from high-throughput sequencing data of multiple 
diploid individuals from a population by a maximum-likelihood (ML) method.
Random union of gametes is assumed in the estimation.  In addition,
all sequence reads are assumed to independently cover at most one
of the two polymorphic sites of interest.

*/

//TODO some performance profiling on buffer sizes.
#define BUFFER_SIZE 1000

using namespace std;

/****************************** estimate_D() *******************************/
/** Function for estimating the LD coefficient D                          **/
/** Returns the ML estimate of D, its minimum and maximum possivle values **/
/**                                                                       **/
/** Parameters:                                                           **/
/** mlNuc1_1: identity of the major allele (1:A, 2:C, 3:G, 4:T)           **/
/** at the first polymorphic site                                         **/
/** mlNuc2_1: identity of the minor allele (1:A, 2:C, 3:G, 4:T)           **/
/** at the first polymorphic site                                         **/
/** mlNuc1_2: identity of the major allele (1:A, 2:C, 3:G, 4:T)           **/
/** at the second polymorphic site                                        **/
/** mlNuc2_2: identity of the minor allele (1:A, 2:C, 3:G, 4:T)           **/
/** at the second polymorphic site                                        **/
/** best_p: ML estimate of the major-allele frequency at the              **/
/** first polymorphic site                                                **/
/** best_q: ML estimate of the major-allele frequency at the              **/
/** second polymorphic site                                               **/                    
/** best_error_1: ML estimate of the error rate at the first              **/
/** polymorphic site                                                      **/
/** best_error_2: ML estimate of the error rate at the second             **/
/** polymorphic site						          **/
/** nsample: number of sampled individuals                                **/
/** mononuc_count_1: nucleotide read quartets at the first                **/
/** polymorphic site                                                      **/
/** mononuc_count_2: nucleotide read quartets at the second               **/
/** polymorphic site                                                      **/
/** cov1: depth of coverage at the first polymorphic site                 **/
/** cov2: depth of coverage at the second polymorphic site                **/
/**                                                                       **/
/** Returns:                                                              **/
/** ML estimates of LD measures                                           **/
/***************************************************************************/
count_t
count_sites(const Locus &X, const Locus &Y)
{
	count_t ret=0;
	std::vector<quartet_t>::const_iterator it_X=X.cbegin();
	std::vector<quartet_t>::const_iterator it_Y=Y.cbegin();
	std::vector<quartet_t>::const_iterator end=Y.cend();
	while (it_Y!=end){
		ret+=(count(*it_X)!=0) & (count(*it_Y)!=0);
		++it_X;
		++it_Y;
	}
	return ret;
}

Linkage estimate_D (const float_t &Ni, const uint8_t &mlNuc1_1, const uint8_t &mlNuc2_1, const uint8_t &mlNuc1_2, const uint8_t &mlNuc2_2, float_t &best_p, float_t &best_q, float_t &best_error_1, float_t &best_error_2, const count_t &nsample, const Locus &mononuc_count_1, const Locus &mononuc_count_2)
{
	Linkage est;
	uint8_t mlNuc3_1, mlNuc4_1, mlNuc3_2, mlNuc4_2;
	float_t maxll, third_best_error_1, prob_mononuc1[11][5], third_best_error_2, prob_mononuc2[11][5];
	float_t size_grid_D, mlD, prob_geno[11];
	int max_mdg, mdg, mgg;
	float_t llhood, null_llhood;
	float_t *prob_obs_mononuc1, *prob_obs_mononuc2, *prob_all_obs;

	prob_obs_mononuc1=new float_t[nsample+1];
	prob_obs_mononuc2=new float_t[nsample+1];
	prob_all_obs=new float_t[nsample+1];
	
	int mig, **ml_mc_1, **ml_mc_2;

	ml_mc_1=new int* [nsample+1];
	for (int i=0; i<nsample+1;++i) 
		ml_mc_1[i]=new int[5];
	ml_mc_2=new int* [nsample+1];
	for (int i=0; i<nsample+1;++i) 
		ml_mc_2[i]=new int[5];

	float_t t_best_D;	// temporarily stores ML estimate of D
	float_t Dmin, Dmax, adj_Dmin, adj_Dmax;

	// Estimate the LD coefficient D between the polymorphic sites
	est.set_Ni(Ni);
	est.set_p(best_p);
	est.set_q(best_q);

	maxll = -FLT_MAX;

	// Find the minimum and maximum possible values of D given the estimated allele frequencies
	if ( best_p*best_q <= (1.0-best_p)*(1.0-best_q) ) {
		Dmin = -best_p*best_q;
	} else {
		Dmin = -(1.0-best_p)*(1.0-best_q);
	}
	if ( best_p*(1.0-best_q) <= (1.0-best_p)*best_q ) {
		Dmax = best_p*(1.0-best_q);
	} else {
		Dmax = (1.0-best_p)*best_q;
	}


	// TODO de-obfuscate. What dose mlNuc1_1 do? 
	// Obtain the nucleotide read counts for the ML analysis
	if (mlNuc1_1*mlNuc2_1 == 2) {
		mlNuc3_1 = 3;
		mlNuc4_1 = 4;
	} else if (mlNuc1_1*mlNuc2_1 == 3) {
		mlNuc3_1 = 2;
		mlNuc4_1 = 4;
	} else if (mlNuc1_1*mlNuc2_1 == 4) {
		mlNuc3_1 = 2;
		mlNuc4_1 = 3;
	} else if (mlNuc1_1*mlNuc2_1 == 6) {
		mlNuc3_1 = 1;
		mlNuc4_1 = 4;
	} else if (mlNuc1_1*mlNuc2_1 == 8) {
		mlNuc3_1 = 1;
		mlNuc4_1 = 3;
	} else if (mlNuc1_1*mlNuc2_1 == 12) {
		mlNuc3_1 = 1;
		mlNuc4_1 = 2;
	}

	if (mlNuc1_2*mlNuc2_2 == 2) {
		mlNuc3_2 = 3;
		mlNuc4_2 = 4;
	} else if (mlNuc1_2*mlNuc2_2 == 3) {
		mlNuc3_2 = 2;
		mlNuc4_2 = 4;
	} else if (mlNuc1_2*mlNuc2_2 == 4) {
		mlNuc3_2 = 2;
		mlNuc4_2 = 3;
	} else if (mlNuc1_2*mlNuc2_2 == 6) {
		mlNuc3_2 = 1;
		mlNuc4_2 = 4;
	} else if (mlNuc1_2*mlNuc2_2 == 8) {
		mlNuc3_2 = 1;
		mlNuc4_2 = 3;
	} else if (mlNuc1_2*mlNuc2_2 == 12) {
		mlNuc3_2 = 1;
		mlNuc4_2 = 2;
	}
	// Calculate the probability for use in the likelihood function
	// Calculate the probability of a nucleotide at the first site given a genotype of the individual
	third_best_error_1 = best_error_1/3.0;
	prob_mononuc1[1][1] = 1.0-best_error_1;			 // A|AB/AB
	prob_mononuc1[1][2] = third_best_error_1;		 // a|AB/AB	
	prob_mononuc1[1][3] = third_best_error_1;		// e1|AB/AB
	prob_mononuc1[1][4] = third_best_error_1;		// e2|AB/AB
	prob_mononuc1[2][1] = 1.0-best_error_1;			//  A|Ab/Ab
	prob_mononuc1[2][2] = third_best_error_1;		//  a|Ab/Ab
	prob_mononuc1[2][3] = third_best_error_1;		// e1|Ab/Ab
	prob_mononuc1[2][4] = third_best_error_1;		// e2|Ab/Ab
	prob_mononuc1[3][1] = third_best_error_1;		//  A|aB/aB
	prob_mononuc1[3][2] = 1.0-best_error_1;			//  a|aB/aB
	prob_mononuc1[3][3] = third_best_error_1;		// e1|aB/aB
	prob_mononuc1[3][4] = third_best_error_1;		// e2|aB/aB
	prob_mononuc1[4][1] = third_best_error_1;		//  A|aB/aB
	prob_mononuc1[4][2] = 1.0-best_error_1;			//  a|aB/aB
	prob_mononuc1[4][3] = third_best_error_1;		// e1|aB/aB
	prob_mononuc1[4][4] = third_best_error_1;		// e2|aB/aB
	prob_mononuc1[5][1] = 1.0-best_error_1;			//  A|AB/Ab
	prob_mononuc1[5][2] = third_best_error_1;		//  a|AB/Ab
	prob_mononuc1[5][3] = third_best_error_1;		// e1|AB/Ab
	prob_mononuc1[5][4] = third_best_error_1;		// e2|AB/Ab
	prob_mononuc1[6][1] = third_best_error_1;		//  A|aB/ab
	prob_mononuc1[6][2] = 1.0-best_error_1;			//  a|aB/ab
	prob_mononuc1[6][3] = third_best_error_1;		// e1|aB/ab
	prob_mononuc1[6][4] = third_best_error_1;		// e2|aB/ab
	prob_mononuc1[7][1] = 0.5 - third_best_error_1/2.;	//  A|AB/aB
	prob_mononuc1[7][2] = 0.5 - third_best_error_1/2.; 	//  a|AB/aB
	prob_mononuc1[7][3] = third_best_error_1;		// e1|AB/aB
	prob_mononuc1[7][4] = third_best_error_1;		// e2|AB/aB
	prob_mononuc1[8][1] = 0.5 - third_best_error_1/2.; 	//  A|Ab/ab
	prob_mononuc1[8][2] = 0.5 - third_best_error_1/2.; 	//  a|Ab/ab
	prob_mononuc1[8][3] = third_best_error_1;		// e1|Ab/ab
	prob_mononuc1[8][4] = third_best_error_1;		// e2|Ab/ab
	prob_mononuc1[9][1] = 0.5 - third_best_error_1/2.;	//  A|AB/ab
	prob_mononuc1[9][2] = 0.5 - third_best_error_1/2.;	//  a|AB/ab
	prob_mononuc1[9][3] = third_best_error_1;		// e1|AB/ab
	prob_mononuc1[9][4] = third_best_error_1;		// e2|AB/ab
	prob_mononuc1[10][1] = 0.5 - third_best_error_1/2.;	//  A|Ab/aB
	prob_mononuc1[10][2] = 0.5 - third_best_error_1/2.;	//  a|Ab/aB
	prob_mononuc1[10][3] = third_best_error_1;		// e1|Ab/aB
	prob_mononuc1[10][4] = third_best_error_1;		// e2|Ab/aB

	// Calculate the probability of a nucleotide at the second site given a genotype of the individual
	third_best_error_2 = best_error_2/3.0;
	prob_mononuc2[1][1] = 1.0-best_error_2;			//  B|AB/AB
	prob_mononuc2[1][2] = third_best_error_2;		//  b|AB/AB
	prob_mononuc2[1][3] = third_best_error_2;		// e1|AB/AB
	prob_mononuc2[1][4] = third_best_error_2;		// e2|AB/AB
	prob_mononuc2[2][1] = third_best_error_2;		//  B|Ab/Ab
	prob_mononuc2[2][2] = 1.0-best_error_2;			//  b|Ab/Ab 
	prob_mononuc2[2][3] = third_best_error_2;		// e1|Ab/Ab
	prob_mononuc2[2][4] = third_best_error_2;		// e2|Ab/Ab
	prob_mononuc2[3][1] = 1.0-best_error_2;			//  B|aB/aB
	prob_mononuc2[3][2] = third_best_error_2;		//  b|aB/aB
	prob_mononuc2[3][3] = third_best_error_2;		// e1|aB/aB
	prob_mononuc2[3][4] = third_best_error_2;		// e2|aB/aB
	prob_mononuc2[4][1] = third_best_error_2;		//  B|ab/ab
	prob_mononuc2[4][2] = 1.0-best_error_2;			//  b|ab/ab
	prob_mononuc2[4][3] = third_best_error_2;		// e1|ab/ab
	prob_mononuc2[4][4] = third_best_error_2;		// e2|ab/ab
	prob_mononuc2[5][1] = 0.5 - third_best_error_2/2.;	//  B|AB/Ab
	prob_mononuc2[5][2] = 0.5 - third_best_error_2/2.;	//  b|AB/Ab
	prob_mononuc2[5][3] = third_best_error_2;		// e1|AB/Ab
	prob_mononuc2[5][4] = third_best_error_2;		// e2|AB/Ab
	prob_mononuc2[6][1] = 0.5 - third_best_error_2/2.;	//  B|aB/ab
	prob_mononuc2[6][2] = 0.5 - third_best_error_2/2.;	//  b|aB/ab
	prob_mononuc2[6][3] = third_best_error_2;		// e1|aB/ab
	prob_mononuc2[6][4] = third_best_error_2;		// e2|aB/ab
	prob_mononuc2[7][1] = 1.0-best_error_2;			//  B|AB/aB
	prob_mononuc2[7][2] = third_best_error_2;		//  b|AB/aB
	prob_mononuc2[7][3] = third_best_error_2;		// e1|AB/aB
	prob_mononuc2[7][4] = third_best_error_2;		// e2|AB/aB
	prob_mononuc2[8][1] = third_best_error_2;		//  B|Ab/ab
	prob_mononuc2[8][2] = 1.0-best_error_2;			//  b|Ab/ab
	prob_mononuc2[8][3] = third_best_error_2;		// e1|Ab/ab
	prob_mononuc2[8][4] = third_best_error_2;		// e2|Ab/ab
	prob_mononuc2[9][1] = 0.5 - third_best_error_2/2.;	//  B|AB/ab
	prob_mononuc2[9][2] = 0.5 - third_best_error_2/2.;	//  b|AB/ab
	prob_mononuc2[9][3] = third_best_error_2;		// e1|AB/ab
	prob_mononuc2[9][4] = third_best_error_2;		// e2|AB/ab
	prob_mononuc2[10][1] = 0.5 - third_best_error_2/2.;	//  B|Ab/aB	
	prob_mononuc2[10][2] = 0.5 - third_best_error_2/2.;	//  b|Ab/aB
	prob_mononuc2[10][3] = third_best_error_2;		// e1|Ab/aB
	prob_mononuc2[10][4] = third_best_error_2;		// e2|Ab/aB

	// Loop over the candidate LD coefficients D between the sites
	size_grid_D = 1.0/(2.0*nsample);

	if ( 2.0*nsample*( Dmax-Dmin ) - (int)( 2.0*nsample*( Dmax-Dmin ) ) >= 0.5) {
		max_mdg = (int)( 2.0*nsample*( Dmax-Dmin ) ) + 2;
	} else {
		max_mdg = (int)( 2.0*nsample*( Dmax-Dmin ) ) + 1;
	}

	for (mdg = 1; mdg <= max_mdg+1; mdg++) {
		if (mdg == max_mdg+1) {
			mlD = 0.0;
		} else {
			mlD = Dmin + (mdg-1)*size_grid_D;
		}
		// Calculate the probability of each of the ten genotypes for an individual
		prob_geno[1] = pow(best_p*best_q+mlD,2.0);					// AB/AB
		prob_geno[2] = pow(best_p*(1.0-best_q)-mlD, 2.0);				// Ab/Ab
		prob_geno[3] = pow((1.0-best_p)*best_q-mlD, 2.0);				// aB/aB
		prob_geno[4] = pow((1.0-best_p)*(1.0-best_q)+mlD, 2.0);				// ab/ab
		prob_geno[5] = 2.0*(best_p*best_q+mlD)*(best_p*(1.0-best_q)-mlD);		// AB/Ab
		prob_geno[6] = 2.0*((1.0-best_p)*best_q-mlD)*((1.0-best_p)*(1.0-best_q)+mlD);	// aB/ab
		prob_geno[7] = 2.0*(best_p*best_q+mlD)*((1.0-best_p)*best_q-mlD);		// AB/aB
		prob_geno[8] = 2.0*(best_p*(1.0-best_q)-mlD)*((1.0-best_p)*(1.0-best_q)+mlD);	// Ab/ab
		prob_geno[9] = 2.0*(best_p*best_q+mlD)*((1.0-best_p)*(1.0-best_q)+mlD);		// AB/ab
		prob_geno[10] = 2.0*(best_p*(1.0-best_q)-mlD)*((1.0-best_p)*best_q-mlD);	// Ab/aB
				
		// Sum the log-likelihoods over the individuals
		llhood = 0.0;

		for (mig = 1; mig <= nsample; mig++) {
			ml_mc_1[mig][1] = mononuc_count_1.get_quartet(mig-1)[mlNuc1_1-1];
			ml_mc_1[mig][2] = mononuc_count_1.get_quartet(mig-1)[mlNuc2_1-1];
			ml_mc_1[mig][3] = mononuc_count_1.get_quartet(mig-1)[mlNuc3_1-1];
			ml_mc_1[mig][4] = mononuc_count_1.get_quartet(mig-1)[mlNuc4_1-1];
			ml_mc_2[mig][1] = mononuc_count_2.get_quartet(mig-1)[mlNuc1_2-1];
			ml_mc_2[mig][2] = mononuc_count_2.get_quartet(mig-1)[mlNuc2_2-1];
			ml_mc_2[mig][3] = mononuc_count_2.get_quartet(mig-1)[mlNuc3_2-1];
			ml_mc_2[mig][4] = mononuc_count_2.get_quartet(mig-1)[mlNuc4_2-1];
			// Sum the probabilities over the genotypes of the individual
			prob_obs_mononuc1[mig] = 0.0;
			prob_obs_mononuc2[mig] = 0.0;
			prob_all_obs[mig] = 0.0;
			count_t cov1=mononuc_count_1.getcoverage(mig-1);
			count_t cov2=mononuc_count_2.getcoverage(mig-1);

			if ( cov1!=0 && cov2!=0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					
					prob_obs_mononuc1[mig] = pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
					prob_obs_mononuc2[mig] = pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
					prob_all_obs[mig] += prob_geno[mgg]*prob_obs_mononuc1[mig]*prob_obs_mononuc2[mig];
				}
				if (prob_all_obs[mig] > 0) {
					llhood = llhood + log(prob_all_obs[mig]);
				} else {
					llhood = -FLT_MAX;
					break;
				}
			} else if (cov1 > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] += pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
				}
				if (prob_obs_mononuc1[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc1[mig]);
				} else {
					llhood = -FLT_MAX;
					break;
				}
			} else if (cov2 > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc2[mig] +=pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
				}
				if (prob_obs_mononuc2[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc2[mig]);
				} else {
					llhood = -FLT_MAX;
					break;
				} 
			} 
		} // End the loop over the individuals
		
		if (mdg == max_mdg+1) {
			null_llhood = llhood;
		} else {
			// Examine whether this is a new ML solution for the sample
			if (llhood > maxll) {
				maxll = llhood;
				t_best_D = mlD;
			}
		}
	}
	
	// OH HELL NO!!!! 
	// if (null_llhood >= maxll) {
	//	maxll = null_llhood;
	// }
	est.set_abs_pos(mononuc_count_1.get_abs_pos() );
	est.set_abs_pos_y(mononuc_count_2.get_abs_pos() );
	est.set_D(t_best_D);
	est.set_fit(maxll);
	est.set_null(null_llhood);

	/*These might need to go global*/
	delete prob_obs_mononuc1;
	delete prob_obs_mononuc2;
	delete prob_all_obs;
	
	for (int i=0; i<nsample+1;++i) 
		delete ml_mc_1[i];
	for (int i=0; i<nsample+1;++i) 
		delete ml_mc_2[i];
	delete ml_mc_1;
	delete ml_mc_2;

	return(est);
}
	
int PopLD(int argc, char *argv[])
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
	env.set_author("Takahiro Maruki");
	env.set_description("Uses a maximum likelihood approach to estimate gametic phase disequilibrium from population data.");

	env.required_arg('p',"pro", 	pro_name,	"please enter a string.", "the input 'pro' file.");
	env.required_arg('m',"map", 	map_name,	"please enter a string.", "the input 'map' file.");
	env.optional_arg('M',"min_n", 	min_number, 	"please enter a number.", "the minimum number of individuals at a site need to calculate LD (default: 10).");
	env.optional_arg('D',"max_d", 	max_d,		"please enter a number.", "the maximum distance between sites for LD calculation (default: none).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.
	
	// Default values of the options

	int argz = 1;	// argument counter

	char *out_file_name;

	Indexed_file <Locus> pro_in;
	Indexed_file <Allele> map_in;

	map_in.open(map_name.c_str() ,std::ios::in);	// Try to open the input file
	pro_in.open(pro_name.c_str(), std::ios::in);	// Try to open the input file
	
	Indexed_file <Linkage> ld_out;	
	ld_out.open(std::ios::out);

	Allele allele1=map_in.read_header();
	Allele allele2=allele1;	
	Locus locus1=pro_in.read_header();	
	Locus locus2=locus1;	
	File_index index=pro_in.get_index();
	ld_out.set_index(index);
	Linkage linkage;
	ld_out.write_header(linkage);

	std::list <Locus> locus_list;

	Locus locus_buffer1[BUFFER_SIZE];
	Locus locus_buffer2[BUFFER_SIZE];

	std::list <Allele> allele_list;

	Allele allele_buffer1[BUFFER_SIZE];
	Allele allele_buffer2[BUFFER_SIZE];

	std::fill_n(locus_buffer1, BUFFER_SIZE, locus1);
	std::fill_n(locus_buffer2, BUFFER_SIZE, locus1);
	std::fill_n(allele_buffer1, BUFFER_SIZE, allele1);
	std::fill_n(allele_buffer2, BUFFER_SIZE, allele1);

	Linkage linkage_buffer[BUFFER_SIZE]; 

	std::fill_n(linkage_buffer, BUFFER_SIZE, linkage);

	count_t nsample=locus1.get_sample_names().size();

	id1_t loc=0, read=0;

	locus_list.insert(locus_list.end(), BUFFER_SIZE, locus1);
	allele_list.insert(allele_list.end(), BUFFER_SIZE, allele1);

	std::list <Locus>::iterator s_locus=locus_list.begin(), e_locus=locus_list.begin(), end_locus=locus_list.end();
	std::list <Allele>::iterator s_allele=allele_list.begin(), e_allele=allele_list.begin(), end_allele=allele_list.end();

	while(e_allele!=end_allele){
		pro_in.read(*(e_locus) );
		map_in.read(*(e_allele) );
		id1_t map_pos=map_in.get_pos(*e_allele);
		while(pro_in.get_pos(*e_locus)<map_pos && pro_in.table_is_open() ) 
				pro_in.read(*e_locus);
		e_locus++;
		e_allele++;
	}

	e_locus=locus_list.begin();
	e_allele=allele_list.begin();
						
	do {
		read=0;
		do {
			size_t number;
			id1_t pos1=locus1.get_abs_pos();
			id1_t pos2=locus2.get_abs_pos();
			id0_t scf1=index.get_id0(pos1);
			id0_t scf2=index.get_id0(pos2);

			//TODO Fix my terrible code! This code block checks to see if the two sites (locus1 and locus2) pass critera for calculating LD.
			if ( (pos2-pos1)<max_d && scf1==scf2) {
				if ( (pos2-pos1)>min_dist ) {
					number=count_sites(locus1, locus2);
					//std::cerr << number << ", " << pos1 << ", " << pos2 << ", " << pos2-pos1 << " ding " << std::distance(e_allele,end_allele) << " -- " << std::distance(s_allele, end_allele) <<  ".\n";
					if ( number >= min_number ) {
						locus_buffer1[read]=locus1;
						allele_buffer1[read]=allele1;
						locus_buffer2[read]=locus2;
						allele_buffer2[read]=allele2;
						read++;
						//size_t Ni=count_sites(locus_buffer1[x], locus_buffer2[x]);
					}
				}
			} else { 
				//If the distance between the two sites is too great, move the first site up.
				s_locus++;
				s_allele++;
				if (s_locus!=end_locus && e_locus!=end_locus){
					locus1=*s_locus;
					allele1=*s_allele;
					e_locus=s_locus;
					e_allele=s_allele;
					locus_list.pop_front();
					allele_list.pop_front();
				} else {
					e_allele=end_allele;
				}
			}

			//Make sure site two isn't at the end of our buffer.
			if (e_allele!=end_allele) {
				//If it isn't, just move site two up..
				locus2=*(e_locus);
				allele2=*(e_allele);
				e_locus++;
				e_allele++;
			} else if (map_in.table_is_open() ) {
				//If it is we have to increase the size of our buffer.
				//std::cerr << "Buffer read\n";
				std::list <Locus>::iterator new_locus=e_locus;
				std::list <Allele>::iterator new_allele=e_allele;
				new_locus--;
				new_allele--;
				locus_list.insert(e_locus, BUFFER_SIZE, locus1);
				allele_list.insert(e_allele, BUFFER_SIZE, allele1);

				end_allele=allele_list.end();
				end_locus=locus_list.end();

				e_locus=new_locus;
				e_allele=new_allele;

				while(new_locus!=end_locus){
					pro_in.read( *(new_locus) );
						map_in.read( *(new_allele) );
					id1_t map_pos=map_in.get_pos(*new_allele);
					while(pro_in.get_pos(*new_locus)<map_pos && pro_in.table_is_open() ) 
						pro_in.read(*new_locus);
					new_locus++;
					new_allele++;
				}

				++e_locus;
				++e_allele;

				locus2=*(e_locus);
				allele2=*(e_allele);
			} else {
				//std::cerr << "Moving locus1\n";
				s_locus++;
				s_allele++;
				if (s_locus!=end_locus){
					locus1=*s_locus;
					allele1=*s_allele;
					e_locus=s_locus;
					e_allele=s_allele;
					locus_list.pop_front();
					allele_list.pop_front();
				} else {
					e_allele=end_allele;
				}
			}
		} while (read<BUFFER_SIZE && (map_in.table_is_open() || s_allele!=end_allele) );
	
		#ifndef NOOMP
		#pragma omp for
		#endif
		for (uint32_t x=0; x<BUFFER_SIZE; ++x){
			if (x<read){
				size_t Ni=count_sites(locus_buffer1[x], locus_buffer2[x]);
			// Estimate the LD coefficient D between the polymorphic sites 
				linkage_buffer[x] = estimate_D(Ni, (uint8_t)allele_buffer1[x].major+1, (uint8_t)allele_buffer1[x].minor+1, (uint8_t)allele_buffer2[x].major+1, (uint8_t)allele_buffer2[x].minor+1, allele_buffer1[x].freq, allele_buffer2[x].freq, allele_buffer1[x].error, allele_buffer2[x].error, nsample, locus_buffer1[x], locus_buffer2[x] );
			}
		}
		for (size_t c=0; c<read; ++c){ 
			ld_out.write(linkage_buffer[c]);
		}
	} while (map_in.table_is_open() || s_allele!=end_allele);
	ld_out.close();
	return 0;
}
