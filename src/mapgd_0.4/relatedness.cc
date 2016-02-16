#include "relatedness.h"

#define BUFFER_SIZE 500

//Moved to inmemory
std::map <Genotype_pair_tuple, size_t> hash_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y)
{
	std::stringstream fb_copy(file_buffer.str() );
	
	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	Population genotypes=gcf_in.read_header();
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		if (genotypes.likelihoods[x].N>1 && genotypes.likelihoods[y].N>1 ){
			counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 2)]+=1;
		}
	}
	return counts;
}

/*Does a regression of allele frequency of the samples on the popualtion allele frequency*/
void set_e(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
	float e;
	e=0;
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
/*	relatedness.f_X_=0;
	relatedness.f_Y_=0;
	relatedness.theta_XY_=0;
	relatedness.gamma_XY_=0;
	relatedness.gamma_YX_=0;

	relatedness.Delta_XY_=0;
	relatedness.delta_XY_=0;*/
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

float_t
get_ll (const Relatedness &rel, const Genotype_pair &pair, const float_t count) 
{
	/*This is the basic likelihood model that we are fitting. It essentially calculates the correlation coefficents
	for the first four moments of the joint distribution. The math needs to be cleaned up here. For now it is typed
	up to minimize the chance of typos. a, b, c and d are r.v. representing the presence (1) or absence of (0) of
	the Major allele in the haploid genomes of idividuals A (for a and b) and C (c and d). A is the r.v. defined as
	A~(a+b)/2 and C~(c+d)/2. Generally what we are doing here is calculating the expectations (e.g. E_A2) of the 
	joint distributions and using that to calculate the joint distribution itself. Problems can occur when 
	correlation coefficents are less than zero, becuase the probabilites of various observations (e.g. mm1mm2) can 
	become negative. These probabilities are forced to be zero, which may be a little arbitrary, but it seems to 
	work.*/


	float_t P, mm1mm2, Mm1mm2, MM1mm2, mm1Mm2, Mm1Mm2, MM1Mm2, mm1MM2, Mm1MM2, MM1MM2;

	P=1-pair.m;

//	std::cerr << P << " is P\n";
	if (P==0) return 0;	
	float_t e=0;
	
	/*This comes from the inverse matrix of the one used to calculate the moments.*/
	
	float_t A=P;//+e(P)*P; 		//mean major allele frequency in A
	float_t C=P;//-e(P)*P; 		//mean major allele frequency in C
	float_t Va=A*(1.-A);	//variance of the two haploid genomes of A 
	float_t Vc=C*(1.-C);	// "   "	"	" 	"     of B
	float_t Sa=sqrt(Va);	//standard deviation of haploid genomes A
	float_t Sc=sqrt(Vc);	// and "	"	"	"	C.

	float_t E_A2  =(rel.f_X_*Va+2.*pow(A,2.)+Va)/2.; //Expectation of A^2
	float_t E_C2  =(rel.f_Y_*Vc+2.*pow(C,2.)+Vc)/2.; //
	float_t E_AC  =rel.theta_XY_*Sa*Sc+A*C;
	float_t ga=(1.-2.*A)/Sa;
	float_t gc=(1.-2.*C)/Sc;
	float_t E_A2C =(rel.gamma_XY_*Va*Sc*ga+A*A*C+Va*(rel.theta_XY_*(1.+2.*A)+rel.f_X_*C+C/(1-C) ) )/2.;
	float_t E_AC2 =(rel.gamma_YX_*Vc*Sa*gc+C*C*A+Vc*(rel.theta_XY_*(1.+2.*C)+rel.f_Y_*A+A/(1-A) ) )/2.;
	float_t ka=1./(1.-A)+1./A-3.;
	float_t kc=1./(1.-C)+1./C-3.;
	float_t E_A2C2=(rel.delta_XY_*sqrt(ka*kc)+rel.Delta_XY_)*Va*Vc+A*A*C*C+rel.f_X_*Va*C*C+rel.f_Y_*Vc*A*A+4.*rel.theta_XY_*Sa*Sc*A*A+C*2.*rel.gamma_XY_*Va*Sc*ga+2.*A*rel.gamma_YX_*Vc*Sa*gc;

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	mm1mm2=1-6*P+0.0*e+2*E_A2+2*E_C2+8.0*E_AC-4*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1mm2=0+4*P+2.0*e-4*E_A2+0*E_C2-10.*E_AC+8*E_A2C+4*E_AC2-2*E_A2C2;
	MM1mm2=0-1*P-0.5*e+2*E_A2+0*E_C2+2.0*E_AC-4*E_A2C-0*E_AC2+1*E_A2C2;
	mm1Mm2=0+4*P-2.0*e+0*E_A2-4*E_C2-10.*E_AC+4*E_A2C+8*E_AC2-2*E_A2C2;
	Mm1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2+12.*E_AC-8*E_A2C-8*E_AC2+4*E_A2C2;
	MM1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+4*E_A2C+0*E_AC2-2*E_A2C2;
	mm1MM2=0-1*P+0.5*e+0*E_A2+2*E_C2+2.0*E_AC+0*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+0*E_A2C+4*E_AC2-2*E_A2C2;
	MM1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2+0.0*E_AC+0*E_A2C+0*E_AC2+1*E_A2C2;


	if (mm1mm2<0 or mm1mm2>1) mm1mm2=0;
	if (Mm1mm2<0 or Mm1mm2>1) Mm1mm2=0; 
	if (MM1mm2<0 or MM1mm2>1) MM1mm2=0; 
	if (mm1Mm2<0 or mm1Mm2>1) mm1Mm2=0; 
	if (Mm1Mm2<0 or Mm1Mm2>1) Mm1Mm2=0;
	if (MM1Mm2<0 or MM1Mm2>1) MM1Mm2=0;
	if (mm1MM2<0 or mm1MM2>1) mm1MM2=0;
	if (Mm1MM2<0 or Mm1MM2>1) Mm1MM2=0;
	if (MM1MM2<0 or MM1MM2>1) MM1MM2=0;

	float_t S=pow(mm1mm2+Mm1mm2+MM1mm2+mm1Mm2+Mm1Mm2+MM1Mm2+mm1MM2+Mm1MM2+MM1MM2, 2);

	if(S>1){
		mm1mm2/=S;
		Mm1mm2/=S;
		MM1mm2/=S;
		mm1Mm2/=S;
		Mm1Mm2/=S;
		MM1Mm2/=S;
		mm1MM2/=S;
		Mm1MM2/=S;
		MM1MM2/=S;
	};
	
	float_t E[9];

        if (mm1mm2>0) E[0]=log(mm1mm2)-pair.X_mm-pair.Y_mm;
     	else E[0]=-FLT_MAX;
	if (Mm1Mm2>0) E[1]=log(Mm1Mm2)-pair.X_Mm-pair.Y_Mm;
     	else E[1]=-FLT_MAX;
        if (MM1MM2>0) E[2]=log(MM1MM2)-pair.X_MM-pair.Y_MM;
     	else E[2]=-FLT_MAX;
        if (mm1Mm2>0) E[3]=log(mm1Mm2)-pair.X_mm-pair.Y_Mm;
     	else E[3]=-FLT_MAX;
        if (MM1Mm2>0) E[4]=log(MM1Mm2)-pair.X_MM-pair.Y_Mm;
     	else E[4]=-FLT_MAX;
        if (Mm1mm2>0) E[5]=log(Mm1mm2)-pair.X_Mm-pair.Y_mm;
     	else E[5]=-FLT_MAX;
        if (Mm1MM2>0) E[6]=log(Mm1MM2)-pair.X_Mm-pair.Y_MM;
     	else E[6]=-FLT_MAX;
        if (mm1MM2>0) E[7]=log(mm1MM2)-pair.X_mm-pair.Y_MM;
     	else E[7]=-FLT_MAX;
        if (MM1mm2>0) E[8]=log(MM1mm2)-pair.X_MM-pair.Y_mm;
     	else E[8]=-FLT_MAX;

	std::sort(E, E+9);

	return (log(1+exp(E[0]-E[8])+exp(E[1]-E[8])+exp(E[2]-E[8])+exp(E[3]-E[8])+exp(E[4]-E[8])+exp(E[5]-E[8])+exp(E[6]-E[8])+exp(E[7]-E[8]) )+E[8] )*count;
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
	
	float_t ll=0;
	
	std::map<Genotype_pair_tuple, size_t>::iterator it=hashed_genotypes_p->begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=hashed_genotypes_p->end();

	while(it!=end){
		first=Genotype_pair::from_tuple(it->first);
		count=it->second;
		ll+=get_ll(rel, first, count);
		if ( isnan(ll) ) break;
		++it;
	}
	if (isnan(ll) ) return FLT_MAX;
	return double(-ll);
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

	gsl_func.f = rel_ll;
	gsl_func.params=&hashed_genotypes;

	s = gsl_multimin_fminimizer_alloc (T, 7);
	gsl_multimin_fminimizer_set (s, &gsl_func, x, ss);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
      
		if (status) break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-5);

		//std::cerr << iter << std::endl;
		/*rel.f_X_ = gsl_vector_get(s->x, 0);
		rel.f_Y_ = gsl_vector_get(s->x, 1);
		rel.theta_XY_ = gsl_vector_get(s->x, 2);
		rel.gamma_XY_ = gsl_vector_get(s->x, 3);
		rel.gamma_YX_ = gsl_vector_get(s->x, 4);
		rel.Delta_XY_ = gsl_vector_get(s->x, 5);
		rel.delta_XY_ = gsl_vector_get(s->x, 6);*/
	}  while (status == GSL_CONTINUE && iter < 600);

	rel.f_X_ = gsl_vector_get(s->x, 0);
	rel.f_Y_ = gsl_vector_get(s->x, 1);
	rel.theta_XY_ = gsl_vector_get(s->x, 2);
	rel.gamma_XY_ = gsl_vector_get(s->x, 3);
	rel.gamma_YX_ = gsl_vector_get(s->x, 4);
	rel.Delta_XY_ = gsl_vector_get(s->x, 5);
	rel.delta_XY_ = gsl_vector_get(s->x, 6);
}


int estimateRel(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

/*	std::string infile="";
	std::vector <size_t> ind;*/
	std::string gcf_name="", rel_name="";

	env_t env;
	env.set_name("mapgd relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	&gcf_name,	&arg_setstr, 	"an error occured while displaying the help message.", "input file name (default stdout)");
	env.optional_arg('o', "output", &rel_name,	&arg_setstr, 	"an error occured while displaying the help message.", "output file name (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;			//The class to read to
	Population genotype;			//The class to write from

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
			relatedness.zero();
			//gestimate(relatedness, hashed_genotypes);
			maximize(relatedness, hashed_genotypes);
			rel_out.write(relatedness);
		}
	}
	rel_out.close();
	return 0;					//Since everything worked, return 0!.
}
