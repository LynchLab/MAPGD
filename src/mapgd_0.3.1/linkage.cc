#include "linkage.h"

/*

Program linkage.cc to estimate linkage disequilibrium (LD) between 
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
};

extern struct estimate estimate_D(int num_pol_sites, int sg, int tg, double N_i, int mlNuc1_1, int mlNuc2_1, int mlNuc1_2, int mlNuc2_2, double best_p, double best_q, double best_error_1, double best_error_2, int nsample, int mononuc_count_1[][5], int mononuc_count_2[][5], int *cov1, int *cov2);

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
struct estimate estimate_D(int num_pol_sites, int sg, int tg, double N_i, int mlNuc1_1, int mlNuc2_1, int mlNuc1_2, int mlNuc2_2, double best_p, double best_q, double best_error_1, double best_error_2, int nsample, int mononuc_count_1[][5], int mononuc_count_2[][5], int *cov1, int *cov2)
{
	vector <estimate> est(num_pol_sites-sg-1);
	int mlNuc3_1, mlNuc4_1, mlNuc3_2, mlNuc4_2;
	double maxll, third_best_error_1, prob_mononuc1[11][5], third_best_error_2, prob_mononuc2[11][5];
	double size_grid_D, mlD, prob_geno[11];
	int max_mdg, mdg, mgg;
	double llhood, prob_obs_mononuc1[nsample+1], prob_obs_mononuc2[nsample+1], prob_all_obs[nsample+1];
	int mig, ml_mc_1[nsample+1][5], ml_mc_2[nsample+1][5];
	double t_best_D;	// temporarily stores ML estimate of D
	double Dmin, Dmax, adj_Dmin, adj_Dmax;

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
	for (mdg = 1; mdg <= max_mdg; mdg++) {
		mlD = Dmin + (mdg-1)*size_grid_D;
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
		
		// Examine whether this is a new ML solution for the sample
		if (llhood > maxll) {
			maxll = llhood;
			t_best_D = mlD;
		}
	
	} // End the loop over the LD coefficients D

	// Calculate the LD measures
	est[tg-sg-1].best_D = t_best_D;
        if (est[tg-sg-1].best_D >= 0) {
        	est[tg-sg-1].best_Dprime = est[tg-sg-1].best_D/Dmax;
        } else {
                est[tg-sg-1].best_Dprime = -est[tg-sg-1].best_D/Dmin;
        }
        est[tg-sg-1].best_D2 = pow(est[tg-sg-1].best_D,2.0);
        est[tg-sg-1].best_r2 = est[tg-sg-1].best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) );
        // Adjust the biases of the LD measures
        est[tg-sg-1].adj_best_D = ( N_i/(N_i-1.0) )*est[tg-sg-1].best_D;
        if (est[tg-sg-1].best_D >= 0) {
        	adj_Dmax = ( N_i/(N_i-1.0) )*Dmax;
                est[tg-sg-1].adj_best_Dprime = est[tg-sg-1].adj_best_D/adj_Dmax;
        } else {
                adj_Dmin = ( N_i/(N_i-1.0) )*Dmin;
                est[tg-sg-1].adj_best_Dprime = -est[tg-sg-1].adj_best_D/adj_Dmin;
        }
	est[tg-sg-1].adj_best_D2 = pow(est[tg-sg-1].adj_best_D,2.0);
	est[tg-sg-1].adj_best_r2 = est[tg-sg-1].adj_best_D2/( best_p*(1.0-best_p)*best_q*(1.0-best_q) ) - 1.0/( (double)N_i );
	// printf("sg: %d\ttg: %d\tbest_D: %f\n", sg, tg, est[tg-sg-1].best_D);
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
	id1_t max_d = INT_MAX;
	my_float_t poly_ll;

	env_t env;
	env.setname("mapgd popld");
	env.setver(VERSION);
	env.setauthor("Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate gametic phase disequalibrium from population data.");

	env.optional_arg('i',"in", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdin).");
	env.optional_arg('o',"out", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file.", "the output file for the program (default stdout).");
	env.optional_arg('m',"max_d",   &max_d,		&arg_setint, 	"an error occured while setting the name of the output file.", "the maximum value fo the distance between polymorphic sites.\n");
	env.optional_arg('l',"polyll",  &poly_ll,	&arg_setfloat_t,"an error occured while setting the name of the output file.", "the minimum value fo the delta ll value of a polymorphic sites.\n");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message.");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version.");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.
	

	if (max_d < INT_MAX) {
		printf("max_d: %d\n", max_d);
	}

	string line;	// String buffer
	profile pro;
	if (infile.size()!=0) {                                 //Iff a filename has been set for infile
		pro.open(infile.c_str(), std::fstream::in);
		if (!pro.is_open() ) {                          //try to open a profile of that name.
			std::cerr << "Cannot open " << infile << " for reading.\n";
			printUsage(env);                        //Print help message on failure.
		}
	}
	else {
		pro.open(std::fstream::in);                     //Iff no filename has been set for infile, open profile from stdin.
	};

	int nsample = pro.size();
	printf("%d individuals to be analyzed\n", nsample);

	// Open the output file and print out the field names
	outstream = fopen(outfile.c_str(), "w");
	fprintf(outstream, "scaffold\tsite 1\tsite 2\tdistance\tbest_D\tbest_D'\tbest_D2\tbest_r2\tadj_best_D\tadj_best_D'\tadj_best_D2\tadj_best_r2\n");
	// printf("scaffold\tsite 1\tsite 2\tdistance\tbest_D\tbest_D'\tbest_D2\tbest_r2\tadj_best_D\tadj_best_D'\tadj_best_D2\tadj_best_r2\n");

	string t_quartet;
	int ig;
	int num_pol_sites;

	std::vector <allele_stat> mle;
	std::vector <Locus> loci;

	while ( pro.read(buffer_site)!=EOF ){
		if(buffer_site.mle.poly_ll>poly_ll){
			loci.push_back(buffer_site.?);
			allele_stat.push_back(buffer_site.mle);
		}
	}

	num_pol_sites = loci.size();
	printf("%d polymorphic sites to be analyzed\n", num_pol_sites);		
	
	int eg, sg;
	string quartet_1[nsample+1];

	for (sg=0; sg<num_pol_sites; sg++) {
	
		vector <estimate> est(num_pol_sites-sg-1);
		est.clear();
		#pragma omp parallel for private(ig, quartet_1) 
		for (int tg=sg+1; tg<num_pol_sites; tg++) {
			#pragma omp critical
			for (ig = 1; ig <= nsample; ig++) {
                        	quartet_1[ig] = quartet[ig].at(sg);
                	}
			string quartet_2[nsample+1];
			int jg, kg, lg, mg;
			int dist_sites, test_dis;
				dist_sites = this_locus->? - that_locus->?;
				if (this_locus->?-that_locus < ?) {
					// printf("dist_sites: %d\n", dist_sites);
					est[tg-sg-1] = estimate_D(num_pol_sites, sg, tg, N_i, mlNuc1.at(sg), mlNuc2.at(sg), mlNuc1.at(tg), mlNuc2.at(tg), best_Maf.at(sg), best_Maf.at(tg), best_error.at(sg), best_error.at(tg), nsample, mononuc_count_1, mononuc_count_2, cov1, cov2);
				}
			}	
		} // End of the loop over the second polymorphic sites
		// Print out the estimation results
		for (eg=0; eg<num_pol_sites-sg-1; eg++) {
                        // Print out the results
			int dist_sites = site.at(eg+sg+1) - site.at(sg);
			if (dist_sites <= max_d) {
				fprintf(outstream, "%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", scaffold.c_str(), site.at(sg), site.at(eg+sg+1), dist_sites, est[eg].best_D, est[eg].best_Dprime, est[eg].best_D2, est[eg].best_r2, est[eg].adj_best_D, est[eg].adj_best_Dprime, est[eg].adj_best_D2, est[eg].adj_best_r2);
			}
		}
	} // End of the loop over the first polymorphic sites
	*/
	return 0;
};
