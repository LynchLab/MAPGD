// Updated on 01/19/16

#include "PopLD.h"

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
	prob_mononuc1[7][1] = 0.5 - third_best_error_1;		//  A|AB/aB
	prob_mononuc1[7][2] = 0.5 - third_best_error_1; 	//  a|AB/aB
	prob_mononuc1[7][3] = third_best_error_1;		// e1|AB/aB
	prob_mononuc1[7][4] = third_best_error_1;		// e2|AB/aB
	prob_mononuc1[8][1] = 0.5 - third_best_error_1; 	//  A|Ab/ab
	prob_mononuc1[8][2] = 0.5 - third_best_error_1; 	//  a|Ab/ab
	prob_mononuc1[8][3] = third_best_error_1;		// e1|Ab/ab
	prob_mononuc1[8][4] = third_best_error_1;		// e2|Ab/ab
	prob_mononuc1[9][1] = 0.5 - third_best_error_1; 	//  A|AB/ab
	prob_mononuc1[9][2] = 0.5 - third_best_error_1;		//  a|AB/ab
	prob_mononuc1[9][3] = third_best_error_1;		// e1|AB/ab
	prob_mononuc1[9][4] = third_best_error_1;		// e2|AB/ab
	prob_mononuc1[10][1] = 0.5 - third_best_error_1;	//  A|Ab/aB
	prob_mononuc1[10][2] = 0.5 - third_best_error_1;	//  a|Ab/aB
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
	prob_mononuc2[5][1] = 0.5 - third_best_error_2;		//  B|AB/Ab
	prob_mononuc2[5][2] = 0.5 - third_best_error_2;		//  b|AB/Ab
	prob_mononuc2[5][3] = third_best_error_2;		// e1|AB/Ab
	prob_mononuc2[5][4] = third_best_error_2;		// e2|AB/Ab
	prob_mononuc2[6][1] = 0.5 - third_best_error_2;		//  B|aB/ab
	prob_mononuc2[6][2] = 0.5 - third_best_error_2;		//  b|aB/ab
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
	prob_mononuc2[9][1] = 0.5 - third_best_error_2;		//  B|AB/ab
	prob_mononuc2[9][2] = 0.5 - third_best_error_2;		//  b|AB/ab
	prob_mononuc2[9][3] = third_best_error_2;		// e1|AB/ab
	prob_mononuc2[9][4] = third_best_error_2;		// e2|AB/ab
	prob_mononuc2[10][1] = 0.5 - third_best_error_2;	//  B|Ab/aB	
	prob_mononuc2[10][2] = 0.5 - third_best_error_2;	//  b|Ab/aB
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
			if (cov1*cov2 > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] = pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
					prob_obs_mononuc2[mig] = pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
					prob_all_obs[mig] = prob_all_obs[mig] + prob_geno[mgg]*prob_obs_mononuc1[mig]*prob_obs_mononuc2[mig];
				}
				if (prob_all_obs[mig] > 0) {
					llhood = llhood + log(prob_all_obs[mig]);
				} else {
					llhood = -FLT_MAX;
				}
			} else if (cov1 > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] = prob_obs_mononuc1[mig] + pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
				}
				if (prob_obs_mononuc1[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc1[mig]);
				} else {
					llhood = -FLT_MAX;
				}
			} else if (cov2 > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc2[mig] = prob_obs_mononuc2[mig] + pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
				}
				if (prob_obs_mononuc2[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc2[mig]);
				} else {
					llhood = -FLT_MAX;
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
	
	// OH HELL NO!!!! You can go to hell. You can go to hell and you can die.
	// if (null_llhood >= maxll) {
	//	maxll = null_llhood;
	// }
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
/*
	std::string infile="";
	std::string outfile="";
	std::string outfilepro;

	bool verbose=false;
	bool quite=false;
	bool noheader=false;
	float_t EMLMIN=0.001;
	count_t MIN=0;
	float_t A=0.00;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	count_t skip=0;
	count_t stop=-1;

	std::vector <int> ind;


	env_t env;
	env.set_name("mapgd ld");
	env.set_version(VERSION);
	env.set_author("Takahiro Maruki");
	env.set_author("Uses a maximum likelihood approach to estimate gametic phase disequalibrium from population data.");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdout).");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file.", "the output file for the program (default stdin).");
	env.optional_arg('p',"out-pro", &outfilepro,	&arg_setstr, 	"an error occured while setting the name of the output file.", "name of a 'cleaned' pro file (default none).");
	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "the individuals to be used in estimates.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).");
	env.optional_arg('m',"minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");

//	env.optional_arg('c',"columns", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "number of columsn in profile (if applicable).");

	env.optional_arg('M',"mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('N',"number", 	&MAXPITCH,	&arg_setint, 	"please provide an int.", "cut-off value for number of bad individuals needed before a site is removed entirely (default 96).");
	env.optional_arg('S',"skip", 	&skip,		&arg_setint, 	"please provide an int.", "number of sites to skip before analysis begins (default 0).");
	env.optional_arg('T',"stop", 	&stop,		&arg_setint, 	"please provide an int.", "maximum number of sites to be analyzed (default All sites)");
	env.flag(	'H',"noheader", &noheader,	&flag_set, 	"takes no argument", "disables printing a headerline.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occured while enabeling verbose excecution.", "prints more information while the command is running.");
	env.flag(	'q', "quite", 	&quite,		&flag_set, 	"an error occured while enabeling quite execution.", "prints less information while the command is running.");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.
	*/
	
	// Default values of the options
	const char* out_file_name = {"Out_PopLD.txt"};
	int max_d = INT_MAX;
	double min_Ni = 10.0;
	int print_help = 0;

	int min_dist=0;
	int argz = 1;	// argument counter

	// Read the specified setting
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-out") == 0) {
			out_file_name = argv[++argz];
		} else if (strcmp(argv[argz], "-max_d") == 0) {
			sscanf(argv[++argz], "%d", &max_d);
		} else if (strcmp(argv[argz], "-min_Ni") == 0) {
                        sscanf(argv[++argz], "%lf", &min_Ni);
		} else {
			fprintf(stderr, "unknown option %s\n", argv[argz]);
			print_help = 1;
			break;
		}
		argz++;
	}
	if (print_help) {	// print error/usage message ?
		fprintf(stderr, "USAGE: %s {<options>}\n", argv[0]);
		fprintf(stderr, "	options: -h: print the usage message\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-max_d <d>: specify the maximum value of the distance between polymorphic sites allowed for estimating LD\n");
		fprintf(stderr, "       -min_Ni <f>: specify the minimum effective sample size required for estimating LD\n");
		exit(1);
	}

	if (max_d != INT_MAX) {
		printf("max_d: %d\n", max_d);
	}

	if (min_Ni != 10.0) {
		printf("min_Ni: %f\n", min_Ni);
	}

	Indexed_file <Allele> in_map1, in_map2;
	Indexed_file <Locus> in_pro1, in_pro2;

	in_map1.open(std::ios::in);	// Try to open the input file
	in_pro1.open(std::ios::in);	// Try to open the input file
	in_map2.open(std::ios::in);	// Try to open the input file
	in_pro2.open(std::ios::in);	// Try to open the input file
	
	Flat_file <Linkage> lsd_out;	// Totally a typo folks

	Allele allele1=in_map1.read_header();
	Allele allele2=allele1;	
	Locus locus1=in_pro1.read_header();	
	Locus locus2=locus1;	

	Linkage linkage;//=

	Locus locus_buffer1[BUFFER_SIZE];
	Allele allele_buffer1[BUFFER_SIZE];

	Locus locus_buffer2[BUFFER_SIZE];
	Allele allele_buffer2[BUFFER_SIZE];

	std::fill_n(locus_buffer1, BUFFER_SIZE, locus1);
	std::fill_n(locus_buffer2, BUFFER_SIZE, locus1);
	std::fill_n(allele_buffer1, BUFFER_SIZE, allele1);
	std::fill_n(allele_buffer2, BUFFER_SIZE, allele1);

	Linkage linkage_buffer[BUFFER_SIZE]; //TODO

	std::fill_n(linkage_buffer, BUFFER_SIZE, linkage);

	count_t nsample=locus1.get_sample_names().size();

	id1_t loc=0, read=0;

//	while ??
	for (size_t x=0; x<BUFFER_SIZE; x++) {
		size_t Ni;
		for (size_t y=x+1; y<BUFFER_SIZE; y++) {
			//locus1=?;
			//locus2=?;
			id1_t pos1=locus1.get_abs_pos();
			id1_t pos2=locus2.get_abs_pos();
			if ( (pos1-pos2)<max_d) {
				if ( (pos1-pos2)>min_dist ) {
					Ni=count_sites(locus_buffer1[x], locus_buffer2[x]);
					if ( Ni >= min_Ni ) {
						locus_buffer1[read]=locus1;
						allele_buffer1[read]=allele1;
						locus_buffer2[read]=locus2;
						allele_buffer2[read]=allele2;
						read++;
					}
				}
			} else break;
		}
	}
	#ifndef NOOMP
	#pragma omp private (loc, read)   
	#endif
			int x=0;
			size_t Ni=count_sites(locus_buffer1[x], locus_buffer2[x]);
			// Estimate the LD coefficient D between the polymorphic sites 
			//TODO dont deleate this
			linkage_buffer[x] = estimate_D(Ni, (uint8_t)allele_buffer1[x].major, (uint8_t)allele_buffer1[x].minor, (uint8_t)allele_buffer2[x].major, (uint8_t)allele_buffer2[x].minor, allele_buffer1[x].freq, allele_buffer2[x].freq, allele_buffer1[x].error, allele_buffer2[x].error, nsample, locus_buffer1[x], locus_buffer2[x] );
//		} 
//	}
	for (size_t c=0; c<read; ++c){ //TODO
		lsd_out.write(linkage_buffer[c]);
	}
	return 0;
}
