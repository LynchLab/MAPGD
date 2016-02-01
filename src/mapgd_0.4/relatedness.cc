#include "relatedness.h"

#define BUFFER_SIZE 500

std::map <Genotype_pair_tuple, size_t> hash_genotypes (Indexed_file <population_genotypes> &gcf_in, size_t x, size_t y)
{
	population_genotypes genotypes;
//	gcf_in.?;
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.is_open() ){
		gcf_in.read(genotypes);
		counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 2)]+=1;
	}
	return counts;
}

/*Does a regression of allele frequency of the samples on the popualtion allele frequency*/
void set_e(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
}

/*Guess starting values of relatedness for the maximization procedure*/
void gestimate(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &counts)
{
	relatedness.zero();
	std::map<Genotype_pair_tuple, size_t>::iterator start=counts.begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=counts.end();
	std::map<Genotype_pair_tuple, size_t>::iterator it=start;
	Genotype_pair pair;

	float_t N=0;
        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
                inc_f(relatedness, pair, it->second);
                inc_theta(relatedness, pair, it->second);
		if (pair.m!=0) N+=it->second/(pair.m*(1-pair.m) );
                ++(it);
        }
	relatedness.f_X_/=N;
	relatedness.f_Y_/=N;
	relatedness.theta_XY_/=N;
	N=0;
	it=start;
        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
               	inc_gamma(relatedness, pair, it->second);
		if (pair.m!=0 && pair.m!=0.5) N+=it->second;
                ++(it);
        }
	relatedness.gamma_XY_/=N;
	relatedness.gamma_YX_/=N;
	it=start;
        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
                inc_Delta(relatedness, pair, it->second);
                ++(it);
        }
	relatedness.Delta_XY_/=N;
}

void inc_f(Relatedness &rel, const Genotype_pair &pair, const size_t &count){
	float_t P=pair.m;
	float_t P2=P*P;
	float_t var=P-P2;
	float_t denom=pow(var, 2);
	if(pair.m!=0){
		rel.f_X_+=(2.*(exp(-pair.X_Mm)/4.+exp(-pair.X_MM)-P2) -var)/denom*count;
		rel.f_Y_+=(2.*(exp(-pair.Y_Mm)/4.+exp(-pair.Y_MM)-P2) -var)/denom*count;
	};
}

void inc_theta(Relatedness &rel, const Genotype_pair &pair, const size_t &count){
	float_t P=pair.m;
	if(P!=0){
		rel.theta_XY_+=( (exp(-pair.X_Mm-pair.Y_Mm)/4.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,2) )/pow(P*(1-P),2)*count;
	};
}

void inc_gamma(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=pair.m;
	if(P!=0 && P!=0.5){
		rel.gamma_XY_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/4.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1-P)*(rel.theta_XY_*(1+2*P)+rel.f_X_*P+P/(1-P) ) )/(P*(1-P) )/(0.5-P)/2*count;
		rel.gamma_YX_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/4.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1.-P)*(rel.theta_XY_*(1.+2.*P)+rel.f_Y_*P+P/(1.-P) ) )/(P*(1.-P) )/(0.5-P)/2.*count;
	};
}

void inc_Delta(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
}

/*Maximizes the relatedness*/
void maximize(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
}


int estimateRel(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

/*	std::string infile="";
	std::vector <size_t> ind;*/
	std::string gcf_name="", rel_name="";

	env_t env;
	env.setname("mapgd relatedness");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	&gcf_name,	&arg_setstr, 	"an error occured while displaying the help message.", "input file name (default stdout)");
	env.optional_arg('o', "output", &rel_name,	&arg_setstr, 	"an error occured while displaying the help message.", "output file name (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <population_genotypes> gcf_in; 	// Open the file with genotypic probabilities.
	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;			//The class to read to
	population_genotypes genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), std::fstream::in);
	else
		gcf_in.open(std::fstream::in);

	if (rel_name.size()!=0)
		rel_out.open(rel_name.c_str(), std::fstream::out);
	else
		rel_out.open(std::fstream::out);

	genotype=gcf_in.read_header();			//This gives us the sample names.
	rel_out.write_header(relatedness);

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;

	for (size_t x=0; x<relatedness.size(); ++x){
		for (size_t y=x+1; x<relatedness.size(); ++x){
			relatedness.set_X_name(genotype.get_sample_names()[x]);
			relatedness.set_Y_name(genotype.get_sample_names()[y]);
			hashed_genotypes=hash_genotypes(gcf_in, x, y);
			set_e(relatedness, hashed_genotypes);
			gestimate(relatedness, hashed_genotypes);
			maximize(relatedness, hashed_genotypes);
			rel_out.write(relatedness);
		}
	}

	return 0;					//Since everything worked, return 0!.
}
