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

using namespace std;

struct estimate
{
	double best_D;		// ML estimate of D
	double best_Dprime;	// ML estimate of D'
	double best_D2;		// ML estimate of D2
	double best_r2;		// ML estimate o r2
	double adj_best_D;	// Bias-adjusted ML estimate of D
	double adj_best_Dprime;	// Bias-adjusted ML estimate of D'
	double adj_best_D2;	// Bias-adjusted ML estimate of D2
	double adj_best_r2;	// Bias-adjusted ML estimate o r2
	double Ni;		// Effective sample size
	double llstat;		// test statistic for examining the statistical significance of LD
};

extern struct estimate estimate_D(int num_pol_sites, int sg, int tg, double Ni, int mlNuc1_1, int mlNuc2_1, int mlNuc1_2, int mlNuc2_2, double best_p, double best_q, double best_error_1, double best_error_2, int nsample, int mononuc_count_1[][5], int mononuc_count_2[][5], int *cov1, int *cov2);

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
struct estimate estimate_D(int num_pol_sites, int sg, int tg, double Ni, int mlNuc1_1, int mlNuc2_1, int mlNuc1_2, int mlNuc2_2, double best_p, double best_q, double best_error_1, double best_error_2, int nsample, int mononuc_count_1[][5], int mononuc_count_2[][5], int *cov1, int *cov2)
{
	vector <estimate> est(num_pol_sites-sg-1);
	int mlNuc3_1, mlNuc4_1, mlNuc3_2, mlNuc4_2;
	double maxll, third_best_error_1, prob_mononuc1[11][5], third_best_error_2, prob_mononuc2[11][5];
	double size_grid_D, mlD, prob_geno[11];
	int max_mdg, mdg, mgg;
	double llhood, prob_obs_mononuc1[nsample+1], prob_obs_mononuc2[nsample+1], prob_all_obs[nsample+1], null_llhood;
	int mig, ml_mc_1[nsample+1][5], ml_mc_2[nsample+1][5];
	double t_best_D;	// temporarily stores ML estimate of D
	double Dmin, Dmax, adj_Dmin, adj_Dmax;

	// printf("Entered the function\n");
	// printf("best_p: %f\tmononuc_count_1[1][1]: %d\tcov1[1]: %d\n", best_p, mononuc_count_1[1][1], cov1[1]);
	// Estimate the LD coefficient D between the polymorphic sites
	est[tg-sg-1].Ni = Ni;
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
			ml_mc_1[mig][1] = mononuc_count_1[mig][mlNuc1_1];
			ml_mc_1[mig][2] = mononuc_count_1[mig][mlNuc2_1];
			ml_mc_1[mig][3] = mononuc_count_1[mig][mlNuc3_1];
			ml_mc_1[mig][4] = mononuc_count_1[mig][mlNuc4_1];
			ml_mc_2[mig][1] = mononuc_count_2[mig][mlNuc1_2];
			ml_mc_2[mig][2] = mononuc_count_2[mig][mlNuc2_2];
			ml_mc_2[mig][3] = mononuc_count_2[mig][mlNuc3_2];
			ml_mc_2[mig][4] = mononuc_count_2[mig][mlNuc4_2];
			// Sum the probabilities over the genotypes of the individual
			prob_obs_mononuc1[mig] = 0.0;
			prob_obs_mononuc2[mig] = 0.0;
			prob_all_obs[mig] = 0.0;
			if (cov1[mig]*cov2[mig] > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] = pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
					prob_obs_mononuc2[mig] = pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
					prob_all_obs[mig] = prob_all_obs[mig] + prob_geno[mgg]*prob_obs_mononuc1[mig]*prob_obs_mononuc2[mig];
				}
				if (prob_all_obs[mig] > 0) {
					llhood = llhood + log(prob_all_obs[mig]);
				} else {
					llhood = -10000000001.0;
				}
			} else if (cov1[mig] > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] = prob_obs_mononuc1[mig] + pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
				}
				if (prob_obs_mononuc1[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc1[mig]);
				} else {
					llhood = -10000000001.0;
				}
			} else if (cov2[mig] > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc2[mig] = prob_obs_mononuc2[mig] + pow(prob_mononuc2[mgg][1],(double)ml_mc_2[mig][1])*pow(prob_mononuc2[mgg][2],(double)ml_mc_2[mig][2])*pow(prob_mononuc2[mgg][3],(double)ml_mc_2[mig][3])*pow(prob_mononuc2[mgg][4],(double)ml_mc_2[mig][4]);
				}
				if (prob_obs_mononuc2[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc2[mig]);
				} else {
					llhood = -10000000001.0;
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
	
	} // End the loop over the LD coefficients D
	if (maxll > -FLT_MAX) {
		// Calculate the LD measures
		est[tg-sg-1].best_D = t_best_D;
        	if (est[tg-sg-1].best_D >= 0) {
        		est[tg-sg-1].best_Dprime = est[tg-sg-1].best_D/Dmax;
        	} else {
                	est[tg-sg-1].best_Dprime = -est[tg-sg-1].best_D/Dmin;
        	}
		if (est[tg-sg-1].best_Dprime > 1.0) {
			est[tg-sg-1].best_Dprime = 1.0;
		} else if (est[tg-sg-1].best_Dprime < -1.0) {
			est[tg-sg-1].best_Dprime = -1.0;
		}
        	est[tg-sg-1].best_D2 = pow(est[tg-sg-1].best_D,2.0);
        	est[tg-sg-1].best_r2 = est[tg-sg-1].best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) );
        	// Adjust the biases of the LD measures
        	est[tg-sg-1].adj_best_D = ( Ni/(Ni-1.0) )*est[tg-sg-1].best_D;
        	if (est[tg-sg-1].best_D >= 0) {
        		adj_Dmax = ( Ni/(Ni-1.0) )*Dmax;
                	est[tg-sg-1].adj_best_Dprime = est[tg-sg-1].adj_best_D/adj_Dmax;
        	} else {
                	adj_Dmin = ( Ni/(Ni-1.0) )*Dmin;
                	est[tg-sg-1].adj_best_Dprime = -est[tg-sg-1].adj_best_D/adj_Dmin;
        	}
		est[tg-sg-1].adj_best_D2 = pow(est[tg-sg-1].adj_best_D,2.0);
		est[tg-sg-1].adj_best_r2 = est[tg-sg-1].adj_best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) ) - 1.0/( (double)Ni );
		// printf("sg: %d\ttg: %d\tbest_D: %f\n", sg, tg, est[tg-sg-1].best_D);
		// Calculate the likelihood-ratio test statistic
		if (null_llhood >= maxll) {
			maxll = null_llhood;
		}
		est[tg-sg-1].llstat = 2.0*(maxll - null_llhood);
	} else {	// ML estimate not found
		est[tg-sg-1].llstat = -FLT_MAX;
	}
	return(est[tg-sg-1]);
}
	
// point to the input and output files
FILE *instream;
FILE *outstream;

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
	env.setname("mapgd ld");
	env.setver(VERSION);
	env.setauthor("Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate gametic phase disequalibrium from population data.");

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

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.
	*/
	
	// Default values of the options
	const char* in_file_name = {"In_PopLD.txt"};
	const char* out_file_name = {"Out_PopLD.txt"};
	int max_d = INT_MAX;
	double min_Ni = 10.0;
	int print_help = 0;

	int argz = 1;	// argument counter

	// Read the specified setting
	while( (argz<argc) && (argv[argz][0] == '-') ) {
		if (strcmp(argv[argz], "-h") == 0) {
			print_help = 1;
		} else if (strcmp(argv[argz], "-in") == 0) {
			in_file_name = argv[++argz];
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
		fprintf(stderr, "	-in <s>: specify the input file name\n");
		fprintf(stderr, "       -out <s>: specify the output file name\n");
		fprintf(stderr, "	-max_d <d>: specify the maximum value of the distance between polymorphic sites allowed for estimating LD\n");
		fprintf(stderr, "       -min_Ni <f>: specify the minimum effective sample size required for estimating LD\n");
		exit(1);
	}

	if (max_d != 2000000000) {
		printf("max_d: %d\n", max_d);
	}

	if (min_Ni != 10.0) {
		printf("min_Ni: %f\n", min_Ni);
	}

	string line;	// String buffer

	ifstream inputFile(in_file_name);	// Try to open the input file
	if ( !inputFile.is_open() ) {	// Exit on failure
		printf("Cannot open %s for reading.\n", in_file_name);
		exit(1);
	}
	
	// Read the header
	string h_scaf, h_site, h_ref_nuc, h_major_allele, h_minor_allele, h_pop_coverage, h_best_p, h_best_error, h_pol_llstat, h_HWE_llstat;	// Stores header names
	getline(inputFile, line);
	istringstream ss(line);
	ss >> h_scaf >> h_site >> h_ref_nuc >> h_major_allele >> h_minor_allele >> h_pop_coverage >> h_best_p >> h_best_error >> h_pol_llstat >> h_HWE_llstat;
	string str;	// Temporarily stores each individual ID
	vector <string> id_ind;      // Stores individual IDs.
	id_ind.clear();
	while (true) {
		ss >> str;
		id_ind.push_back(str);
		if ( ss.eof() ) {
			break;
		}
	}
	int nsample = id_ind.size();
	printf("%d individuals to be analyzed\n", nsample);

	// Open the output file and print out the field names
	outstream = fopen(out_file_name, "w");
	fprintf(outstream, "scaffold\tsite_1\tsite_2\tdistance\tbest_D\tbest_D'\tbest_D2\tbest_r2\tadj_best_D\tadj_best_D'\tadj_best_D2\tadj_best_r2\tNi\tllstat\n");
	// printf("scaffold\tsite_1\tsite_2\tdistance\tbest_D\tbest_D'\tbest_D2\tbest_r2\tadj_best_D\tadj_best_D'\tadj_best_D2\tadj_best_r2\tNi\tllstat\n");

	int t_site;	// temporarily stores the coordinate
	double t_best_Maf, t_best_error;	// temporarily stores major-allele frequency estimates and error rates
	double pol_llstat, HWE_llstat;
	string t_quartet;
	vector <string> quartet[nsample+1];
	vector <int> site;
	vector <double> best_Maf;
	vector <double> best_error;
	int ig;
	string scaffold;
	string ref_nuc;
	string s_mlNuc1, s_mlNuc2;
	int t_mlNuc1, t_mlNuc2;
	int pop_cov;
	vector <int> mlNuc1, mlNuc2;
	int num_pol_sites;

	// Read the main data
	site.clear();
	mlNuc1.clear();
	mlNuc2.clear();
	best_Maf.clear();
	best_error.clear();
	for (ig = 1; ig <= nsample; ig++) {
		quartet[ig].clear();
	}
	while ( getline(inputFile, line) ) {
		istringstream ss(line);
		ss >> scaffold >> t_site >> ref_nuc >> s_mlNuc1 >> s_mlNuc2 >> pop_cov >> t_best_Maf >> t_best_error >> pol_llstat >> HWE_llstat;
		if (s_mlNuc1 == "A") {
			t_mlNuc1 = 1;
		} else if (s_mlNuc1 == "C") {
			t_mlNuc1 = 2;
		} else if (s_mlNuc1 == "G") {
                        t_mlNuc1 = 3;
		} else if (s_mlNuc1 == "T") {
                        t_mlNuc1 = 4;
		} else {
			fprintf(stderr, "Problem in the major-allele designation found at site %d on scaffold %s\n", t_site, scaffold.c_str());
		}
		if (s_mlNuc2 == "A") {
                        t_mlNuc2 = 1;
                } else if (s_mlNuc2 == "C") {
                        t_mlNuc2 = 2;
                } else if (s_mlNuc2 == "G") {
                        t_mlNuc2 = 3;
                } else if (s_mlNuc2 == "T") {
                        t_mlNuc2 = 4;
                } else {
                        fprintf(stderr, "Problem in the minor-allele designation found at site %d on scaffold %s\n", t_site, scaffold.c_str());
                }
		site.push_back(t_site);
		mlNuc1.push_back(t_mlNuc1);
		mlNuc2.push_back(t_mlNuc2);
		best_Maf.push_back(t_best_Maf);
		best_error.push_back(t_best_error);
		for (ig = 1; ig <= nsample; ig++) {
			ss >> t_quartet;
			quartet[ig].push_back(t_quartet);
		}
	}
	num_pol_sites = best_Maf.size();
	printf("%d polymorphic sites to be analyzed\n", num_pol_sites);		
	
	int eg, sg;
	string quartet_1[nsample+1];

	for (sg=0; sg<num_pol_sites; sg++) {
	
		vector <estimate> est(num_pol_sites-sg-1);
		est.clear();
		#ifdef PRAGMA
		#pragma omp parallel for private(ig, quartet_1) 
		#endif
		for (int tg=sg+1; tg<num_pol_sites; tg++) {
			#ifdef PRAGMA
			#pragma omp critical
			#endif
			for (ig = 1; ig <= nsample; ig++) {
                        	quartet_1[ig] = quartet[ig].at(sg);
                	}
			string quartet_2[nsample+1];
			int jg, kg, lg, mg;
			int dist_sites, test_dis, mononuc_count_1[nsample+1][5], mononuc_count_2[nsample+1][5];
			int cov1[nsample+1], cov2[nsample+1];
			double ind_Ni, Ni;
			int size_quartet_1, size_quartet_2, count_nuc, digit, num_digit[5];
			// printf("sg: %d\ttg: %d\n", sg, tg);
			test_dis = 1;
			// printf("Inside the loop\n");
			if (test_dis == 1) {
				for (ig = 1; ig <= nsample; ig++) {
					quartet_2[ig] = quartet[ig].at(tg);
				}
				// printf("site_1: %d\tbest_p: %f\tsite_2: %d\tbest_q: %f\n", site.at(sg), best_Maf.at(sg), site.at(tg), best_Maf.at(tg));
				dist_sites = site.at(tg) - site.at(sg);
				// printf("sg: %d\ttg: %d\tdis_sites: %d\n", sg, tg, dist_sites);
				if (dist_sites > max_d) {
					test_dis = 0;
				}
				if (test_dis == 1) {
					// printf("dist_sites: %d\n", dist_sites);
					Ni = 0.0;
					for (ig = 1; ig <= nsample; ig++) {
						size_quartet_1 = quartet_1[ig].size();
                                        	jg = 1;
                                        	digit = 0;
                                        	for (kg = 0; kg < size_quartet_1; kg++) {
                                                	if ( quartet_1[ig].at(kg) == '/') {
                                                        	mg = pow(10,digit-1);
                                                        	count_nuc = 0;
                                                        	lg = 0;
                                                        	while (lg <= digit-1) {
                                                                	count_nuc = count_nuc + mg*num_digit[lg];
                                                                	mg = mg/10;
                                                                	num_digit[lg] = 0;
                                                                	lg = lg + 1;
                                                        	}
                                                        	mononuc_count_1[ig][jg] = count_nuc;
                                                        	jg = jg + 1;
                                                        	count_nuc = 0;
                                                        	digit = 0;
                                                	} else {
                                                        	num_digit[digit] = quartet_1[ig].at(kg) - '0';
                                                        	digit = digit + 1;
                                                	}
                                                	if (kg == size_quartet_1 - 1) {
                                                        	mg = pow(10,digit-1);
                                                        	count_nuc = 0;
                                                        	lg = 0;
                                                        	while (lg <= digit-1) {
                                                                	count_nuc = count_nuc + mg*num_digit[lg];
                                                                	mg = mg/10;
                                                                	num_digit[lg] = 0;
                                                                	lg = lg + 1;
                                                        	}
                                                        	mononuc_count_1[ig][jg] = count_nuc;
                                                        	jg = 1;
                                                	}
                                        	}
                                        	size_quartet_2 = quartet_2[ig].size();
                                        	jg = 1;
                                        	digit = 0;
                                        	for (kg = 0; kg < size_quartet_2; kg++) {
                                                	if ( quartet_2[ig].at(kg) == '/') {
                                                        	mg = pow(10,digit-1);
                                                        	count_nuc = 0;
                                                        	lg = 0;
                                                        	while (lg <= digit-1) {
                                                                	count_nuc = count_nuc + mg*num_digit[lg];
                                                                	mg = mg/10;
                                                                	num_digit[lg] = 0;
                                                                	lg = lg + 1;
                                                        	}
                                                        	mononuc_count_2[ig][jg] = count_nuc;
                                                        	jg = jg + 1;
                                                        	count_nuc = 0;
                                                        	digit = 0;
                                                	} else {
                                                        	num_digit[digit] = quartet_2[ig].at(kg) - '0';
                                                        	digit = digit + 1;
                                                	}
                                                	if (kg == size_quartet_2 - 1) {
                                                        	mg = pow(10,digit-1);
                                                        	count_nuc = 0;
                                                        	lg = 0;
								while (lg <= digit-1) {
                                                                	count_nuc = count_nuc + mg*num_digit[lg];
                                                                	mg = mg/10;
                                                                	num_digit[lg] = 0;
                                                                	lg = lg + 1;
                                                        	}
                                                        	mononuc_count_2[ig][jg] = count_nuc;
                                                        	jg = 1;
                                                	}
                                        	}
                                        	cov1[ig] = mononuc_count_1[ig][1] + mononuc_count_1[ig][2] + mononuc_count_1[ig][3] + mononuc_count_1[ig][4];
                                        	cov2[ig] = mononuc_count_2[ig][1] + mononuc_count_2[ig][2] + mononuc_count_2[ig][3] + mononuc_count_2[ig][4];
                                        	ind_Ni = ( 1.0-pow( 0.5,(double)cov1[ig] ) )*pow( 0.5,(double)(cov2[ig]+1.0) ) + pow( 0.5,(double)(cov1[ig]+1.0) )*( 1.0-pow(0.5,(double)cov2[ig]) ) + ( 1.0-pow(0.5,(double)cov1[ig]) )*( 1.0-pow(0.5,(double)cov2[ig]) );
                                        	Ni = Ni + ind_Ni;
                                	} // End of the loop over the individuals 
					// printf("Exited from the loop\n");
					if (Ni >= min_Ni) {
						// Estimate the LD coefficient D between the polymorphic sites 
						est[tg-sg-1] = estimate_D(num_pol_sites, sg, tg, Ni, mlNuc1.at(sg), mlNuc2.at(sg), mlNuc1.at(tg), mlNuc2.at(tg), best_Maf.at(sg), best_Maf.at(tg), best_error.at(sg), best_error.at(tg), nsample, mononuc_count_1, mononuc_count_2, cov1, cov2);
					} else {
						est[tg-sg-1].Ni = Ni;
					}
				}
			}	
		} // End of the loop over the second polymorphic sites
		// Print out the estimation results
		for (eg=0; eg<num_pol_sites-sg-1; eg++) {
                        // Print out the results
			int dist_sites = site.at(eg+sg+1) - site.at(sg);
			if (dist_sites <= max_d) {
				if (est[eg].Ni >= min_Ni) {
					if (est[eg].llstat != -FLT_MAX) {
						// printf("%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site.at(sg), site.at(eg+sg+1), dist_sites, est[eg].best_D, est[eg].best_Dprime, est[eg].best_D2, est[eg].best_r2, est[eg].adj_best_D, est[eg].adj_best_Dprime, est[eg].adj_best_D2, est[eg].adj_best_r2, est[eg].Ni, est[eg].llstat);
						fprintf(outstream, "%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site.at(sg), site.at(eg+sg+1), dist_sites, est[eg].best_D, est[eg].best_Dprime, est[eg].best_D2, est[eg].best_r2, est[eg].adj_best_D, est[eg].adj_best_Dprime, est[eg].adj_best_D2, est[eg].adj_best_r2, est[eg].Ni, est[eg].llstat);
					} else {
						fprintf(stderr, "ML estimates not found for sites %d and %d on %s\n", site.at(sg), site.at(eg+sg+1), scaffold.c_str());
						fprintf(outstream, "%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%f\tNA\n", scaffold.c_str(), site.at(sg), site.at(eg+sg+1), dist_sites, est[eg].Ni);
					}
				} else {
					fprintf(outstream, "%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t%f\tNA\n", scaffold.c_str(), site.at(sg), site.at(eg+sg+1), dist_sites, est[eg].Ni);
				}
			}
		}
	} // End of the loop over the first polymorphic sites
	
	return 0;
};
