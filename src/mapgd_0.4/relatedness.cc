#include "relatedness.h"

#define BUFFER_SIZE 500

//Moved to inmemory
std::map <Genotype_pair_tuple, size_t> hash_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y)
{
	std::cerr << "about to reading genotypes\n";
	std::stringstream fb_copy(file_buffer.str() );
	
	Indexed_file <population_genotypes> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	population_genotypes genotypes=gcf_in.read_header();
	std::map <Genotype_pair_tuple, size_t> counts;
	std::cerr << "reading genotypes\n";
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 2)]+=1;
	}
	std::cerr << "done genotypes\n";
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

void 
inc_Delta (Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{

}

double
rel_ll (const gsl_vector *v, void *void_hashed_genotypes_p)
{
	Relatedness rel;
	std::map <Genotype_pair_tuple, size_t> *hashed_genotypes_p=(std::map <Genotype_pair_tuple, size_t> *) void_hashed_genotypes_p;
	Genotype_pair first;
	size_t count;

	rel.f_X_ = gsl_vector_get(v, 0);
	rel.f_Y_ = gsl_vector_get(v, 1);
	rel.theta_XY_ = gsl_vector_get(v, 2);
	rel.gamma_XY_ = gsl_vector_get(v, 3);
	rel.gamma_YX_ = gsl_vector_get(v, 4);
	rel.Delta_XY_ = gsl_vector_get(v, 5);
	rel.delta_XY_ = gsl_vector_get(v, 6);
	
	double ll=0;
	
	std::map<Genotype_pair_tuple, size_t>::iterator it=hashed_genotypes_p->begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=hashed_genotypes_p->end();

	while(it!=end){
		first=Genotype_pair::from_tuple(it->first);
		count=it->second;
		ll+=log((rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) + rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) + pow(first.m, 4) - pow(first.m, 3)*(first.m - 1)*(rel.f_X_ + rel.f_Y_ + 4*rel.theta_XY_) + first.m*(2*rel.gamma_XY_ + 2*rel.gamma_YX_)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m))*exp(-first.X_MM)*exp(-first.Y_MM) + (rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) + rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) - first.m*pow(first.m - 1, 3)*(rel.f_X_ + rel.f_Y_ + 4*rel.theta_XY_) + (2*rel.gamma_XY_ + 2*rel.gamma_YX_)*(first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) + pow(first.m - 1, 4))*exp(-first.X_mm)*exp(-first.Y_mm) + (4*rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) + 4*rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) + 4*pow(first.m, 2)*pow(first.m - 1, 2) - first.m*(first.m - 1)*(4*rel.f_X_*first.m*(first.m - 1) + 4*rel.f_Y_*first.m*(first.m - 1) + 4*rel.theta_XY_*(pow(first.m, 2) + 2*first.m*(first.m - 1) + pow(first.m - 1, 2))) + (4*rel.gamma_XY_ + 4*rel.gamma_YX_)*(2*first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m))*exp(-first.X_Mm)*exp(-first.Y_Mm) - (2*rel.gamma_XY_*first.m*(-2*pow(first.m, 3) + 3*pow(first.m, 2) - first.m) + 2*rel.gamma_YX_*(first.m - 1)*(-2*pow(first.m, 3) + 3*pow(first.m, 2) - first.m) - rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) - pow(first.m, 2)*pow(first.m - 1, 2) + first.m*(first.m - 1)*(rel.f_X_*pow(first.m, 2) + rel.f_Y_*pow(first.m - 1, 2) + 4*rel.theta_XY_*first.m*(first.m - 1)))*exp(-first.Y_MM)*exp(-first.X_mm) + (-4*rel.gamma_XY_*first.m*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.gamma_YX_*(2*first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - 2*rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) + pow(first.m, 3)*(-2*first.m + 2) + first.m*(first.m - 1)*(2*rel.f_X_*pow(first.m, 2) + 2*rel.f_Y_*first.m*(first.m - 1) + 4*rel.theta_XY_*(pow(first.m, 2) + first.m*(first.m - 1))))*exp(-first.Y_MM)*exp(-first.X_Mm) - (2*rel.gamma_XY_*(first.m - 1)*(-2*pow(first.m, 3) + 3*pow(first.m, 2) - first.m) + 2*rel.gamma_YX_*first.m*(-2*pow(first.m, 3) + 3*pow(first.m, 2) - first.m) - rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) - pow(first.m, 2)*pow(first.m - 1, 2) + first.m*(first.m - 1)*(rel.f_X_*pow(first.m - 1, 2) + rel.f_Y_*pow(first.m, 2) + 4*rel.theta_XY_*first.m*(first.m - 1)))*exp(-first.X_MM)*exp(-first.Y_mm) + (-4*rel.gamma_XY_*(first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.gamma_YX_*(2*first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - 2*rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) - 2*first.m*pow(first.m - 1, 3) + first.m*(first.m - 1)*(2*rel.f_X_*pow(first.m - 1, 2) + 2*rel.f_Y_*first.m*(first.m - 1) + 4*rel.theta_XY_*(first.m*(first.m - 1) + pow(first.m - 1, 2))))*exp(-first.X_Mm)*exp(-first.Y_mm) + (-2*rel.gamma_XY_*(2*first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 4*rel.gamma_YX_*first.m*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - 2*rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) + pow(first.m, 3)*(-2*first.m + 2) + first.m*(first.m - 1)*(2*rel.f_X_*first.m*(first.m - 1) + 2*rel.f_Y_*pow(first.m, 2) + 4*rel.theta_XY_*(pow(first.m, 2) + first.m*(first.m - 1))))*exp(-first.X_MM)*exp(-first.Y_Mm) + (-2*rel.gamma_XY_*(2*first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 4*rel.gamma_YX_*(first.m - 1)*(2*pow(first.m, 3) - 3*pow(first.m, 2) + first.m) - 2*rel.delta_XY_*(-6*pow(first.m, 4) + 12*pow(first.m, 3) + 3*pow(first.m, 2)*pow(first.m - 1, 2) - 7*pow(first.m, 2) + first.m) - 2*rel.Delta_XY_*pow(first.m, 2)*pow(first.m - 1, 2) - 2*first.m*pow(first.m - 1, 3) + first.m*(first.m - 1)*(2*rel.f_X_*first.m*(first.m - 1) + 2*rel.f_Y_*pow(first.m - 1, 2) + 4*rel.theta_XY_*(first.m*(first.m - 1) + pow(first.m - 1, 2))))*exp(-first.Y_Mm)*exp(-first.X_mm))
*count;
	
	if ( isnan(ll) ) break;
	++it;
	}
	return ll;
}

/*Maximizes the relatedness*/
void 
maximize(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
	const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function gsl_func;
	
	size_t iter = 0;
	int status;
	double size;

	/* Starting point and stepsizes. I would really prefer to do this with a 
	 * Newton-Raphson method, which just inexplicably like more than the 
	 * Nelder-Mead, but I'm being lasy today, or more accurately there are 
	 * fundemental problems with setting the problem up to use the NR method.
	 */
	x=gsl_vector_alloc(7);
	gsl_vector_set(x, 0, rel.f_X_);
	gsl_vector_set(x, 1, rel.f_Y_);
	gsl_vector_set(x, 2, rel.theta_XY_);
	gsl_vector_set(x, 3, rel.gamma_XY_);
	gsl_vector_set(x, 4, rel.gamma_YX_);
	gsl_vector_set(x, 5, rel.Delta_XY_);
	gsl_vector_set(x, 6, rel.delta_XY_);

	ss=gsl_vector_alloc(7);
	gsl_vector_set_all(ss,0.25);	

	gsl_func.n=7;
	std::cerr << "HI1\n";
	gsl_func.f = rel_ll;
	gsl_func.params=&hashed_genotypes;

	std::cerr << "HI2\n";

	s = gsl_multimin_fminimizer_alloc (T, 7);
	gsl_multimin_fminimizer_set (s, &gsl_func, x, ss);


	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
      
		if (status) break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);

		if (status == GSL_SUCCESS) printf ("converged to minimum at\n");

		printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
			iter,
 			gsl_vector_get (s->x, 0), 
			gsl_vector_get (s->x, 1), 
			s->fval, size);
	}  while (status == GSL_CONTINUE && iter < 100);

	rel.f_X_ = gsl_vector_get(x, 0);
	rel.f_Y_ = gsl_vector_get(x, 1);
	rel.theta_XY_ = gsl_vector_get(x, 2);
	rel.gamma_XY_ = gsl_vector_get(x, 3);
	rel.gamma_YX_ = gsl_vector_get(x, 4);
	rel.Delta_XY_ = gsl_vector_get(x, 5);
	rel.delta_XY_ = gsl_vector_get(x, 6);
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
	Indexed_file <population_genotypes> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;			//The class to read to
	population_genotypes genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), std::ios::in);
	else
		gcf_in.open(std::ios::in);

	if (rel_name.size()!=0)
		rel_out.open(rel_name.c_str(), std::ios::out);
	else
		rel_out.open(std::ios::out);

	genotype=gcf_in.read_header();			//This gives us the sample names.

	gcf_mem.open(&file_buffer, std::ios::out );
	gcf_mem.set_index(gcf_in.get_index() );
	gcf_mem.write_header(genotype);

	while(gcf_in.table_is_open() ){
		gcf_in.read(genotype);
		gcf_mem.write(genotype);
	}
	gcf_mem.close();

	rel_out.write_header(relatedness);

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();

	for (size_t x=0; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
			relatedness.set_X_name(genotype.get_sample_names()[x]);
			relatedness.set_Y_name(genotype.get_sample_names()[y]);
			hashed_genotypes=hash_genotypes(file_buffer, x, y);
			set_e(relatedness, hashed_genotypes);
			gestimate(relatedness, hashed_genotypes);
			maximize(relatedness, hashed_genotypes);
			rel_out.write(relatedness);
		}
	}

	return 0;					//Since everything worked, return 0!.
}
