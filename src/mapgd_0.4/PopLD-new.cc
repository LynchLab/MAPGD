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
	float_t best_D;		// ML estimate of D
	float_t best_Dprime;	// ML estimate of D'
	float_t best_D2;		// ML estimate of D2
	float_t best_r2;		// ML estimate o r2
	float_t adj_best_D;	// Bias-adjusted ML estimate of D
	float_t adj_best_Dprime;	// Bias-adjusted ML estimate of D'
	float_t adj_best_D2;	// Bias-adjusted ML estimate of D2
	float_t adj_best_r2;	// Bias-adjusted ML estimate o r2
	float_t Ni;		// Effective sample size
	float_t llstat;		// test statistic for examining the statistical significance of LD
};

//extern struct estimate estimate_D(int num_pol_sites, int sg, int tg, float_t Ni, int mlNuc1_1, int mlNuc2_1, int mlNuc1_2, int mlNuc2_2, float_t best_p, float_t best_q, float_t best_error_1, float_t best_error_2, int nsample, int mononuc_count_1[][5], int mononuc_count_2[][5], int *cov1, int *cov2);


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
estimate estimate_D (const float_t &Ni, const gt_t &mlNuc1_1, const gt_t &mlNuc2_1, const gt_t &mlNuc1_2, const gt_t &mlNuc2_2, float_t &best_p, float_t &best_q, float_t &best_error_1, float_t &best_error_2, const count_t &nsample, const quartet &mononuc_count_1, const quartet &mononuc_count_2, const count_t *cov1, const count_t *cov2)
{
	estimate est;
	gt_t mlNuc3_1, mlNuc4_1, mlNuc3_2, mlNuc4_2;
	float_t maxll, third_best_error_1, prob_mononuc1[11][5], third_best_error_2, prob_mononuc2[11][5];
	float_t size_grid_D, mlD, prob_geno[11];
	int max_mdg, mdg, mgg;
	float_t llhood, null_llhood;
	float_t *prob_obs_mononuc1, *prob_obs_mononuc2, *prob_all_obs;

	*prob_obs_mononuc1=new float_t[nsample+1];
	*prob_obs_mononuc2=new float_t[nsample+1];
	*prob_all_obs=new float_t[nsample+1];
	
	int mig, ml_mc_1[nsample+1][5], ml_mc_2[nsample+1][5];
	float_t t_best_D;	// temporarily stores ML estimate of D
	float_t Dmin, Dmax, adj_Dmin, adj_Dmax;

	// printf("Entered the function\n");
	// printf("best_p: %f\tmononuc_count_1[1][1]: %d\tcov1[1]: %d\n", best_p, mononuc_count_1[1][1], cov1[1]);
	// Estimate the LD coefficient D between the polymorphic sites
	est.Ni = Ni;
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
					llhood = -FLT_MAX;
				}
			} else if (cov1[mig] > 0) {
				for(mgg = 1; mgg <= 10; mgg++) {
					prob_obs_mononuc1[mig] = prob_obs_mononuc1[mig] + pow(prob_mononuc1[mgg][1],(double)ml_mc_1[mig][1])*pow(prob_mononuc1[mgg][2],(double)ml_mc_1[mig][2])*pow(prob_mononuc1[mgg][3],(double)ml_mc_1[mig][3])*pow(prob_mononuc1[mgg][4],(double)ml_mc_1[mig][4]);
				}
				if (prob_obs_mononuc1[mig] > 0) {
					llhood = llhood + log(prob_obs_mononuc1[mig]);
				} else {
					llhood = -FLT_MAX;
				}
			} else if (cov2[mig] > 0) {
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
	
	} // End the loop over the LD coefficients D
	if (maxll > -FLT_MAX) {
		// Calculate the LD measures
		est.best_D = t_best_D;
        	if (est.best_D >= 0) {
        		est.best_Dprime = est.best_D/Dmax;
        	} else {
                	est.best_Dprime = -est.best_D/Dmin;
        	}
		if (est.best_Dprime > 1.0) {
			est.best_Dprime = 1.0;
		} else if (est.best_Dprime < -1.0) {
			est.best_Dprime = -1.0;
		}
        	est.best_D2 = pow(est.best_D,2.0);
        	est.best_r2 = est.best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) );
        	// Adjust the biases of the LD measures
        	est.adj_best_D = ( Ni/(Ni-1.0) )*est.best_D;
        	if (est.best_D >= 0) {
        		adj_Dmax = ( Ni/(Ni-1.0) )*Dmax;
                	est.adj_best_Dprime = est.adj_best_D/adj_Dmax;
        	} else {
                	adj_Dmin = ( Ni/(Ni-1.0) )*Dmin;
                	est.adj_best_Dprime = -est.adj_best_D/adj_Dmin;
        	}
		est.adj_best_D2 = pow(est.adj_best_D,2.0);
		est.adj_best_r2 = est.adj_best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) ) - 1.0/( (double)Ni );
		// printf("sg: %d\ttg: %d\tbest_D: %f\n", sg, tg, est.best_D);
		// Calculate the likelihood-ratio test statistic
		if (null_llhood >= maxll) {
			maxll = null_llhood;
		}
		est.llstat = 2.0*(maxll - null_llhood);
	} else {	// ML estimate not found
		est.llstat = -FLT_MAX;
	}
	return(est);
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

	if (max_d != INT_MAX) {
		printf("max_d: %d\n", max_d);
	}

	if (min_Ni != 10.0) {
		printf("min_Ni: %f\n", min_Ni);
	}

	Indexed_file <allele_stat> in_map;
	Indexed_file <Locus> in_pro;

	in_map.open(in_file_name.c_str(), std::ios::in);	// Try to open the input file
	in_pro.open(in_file_name.c_str(), std::ios::in);	// Try to open the input file
	
	Indexed_file <Linkage_stat> lds_out;
	
	Locus locus_buffer[BUFFER_SIZE];
	allele_stat allele_buffer[BUFFER_SIZE];

	linkage_data linkage_buffer[?]; //TODO

	num_pol_sites? //TODO

	while (x<BUFFER_SIZE && in_map.is_open() ){
		while (pro_in.get_pos(locus)<in_map.get_pos(allele) ){
				pro_in.read(locus);
		}
		if (pro_in.get_pos(locus)==in_map.get_pos(allele) {
			locus_buffer[x]=locus; 
			allele_buffer[x]=allele; 
		}
		in_map.read(allele);
		++x;
	}
	for (size_t sg=0; sg<BUFFER_SIZE; sg++) {
		
		#ifndef NOOMP
		#pragma omp parallel for private(x, quartet_1) 
		#endif
		for (size_t tg=sg+1; tg<num_pol_sites; tg++) {
			Ni=const_sites(site1, site2);
			dist_sites = pro_in.get_pos(site2) - pro_in.get_pos(site1);
			if (dist_sites <= max_d && Ni >= min_Ni) {
				// Estimate the LD coefficient D between the polymorphic sites 
				//TODO dont deleate this
				linkage_buffer[tg-sg-1] = estimate_D(Ni, allele_buffer[sg].major, allele_buffer[sg].minor, allele_buffer[tg].major, allele_buffer[tg].minor, allele_buffer[sg].freq, allele_buffer[tg].freq, allele_bffuer[sg].error, allele_buffer[tg].error, nsample <-FIXIT, locus_buffer[sg].get_quartets(), locus_buffer[tg].get_quartetes(), locus_buffer[sg].getcount(), locus_buffer[tg].getcount());
			} 
		}
	}	
	for (??){ //TODO
		lds_out.write(linkage_buffer[x]);
	}
	return 0;
};
