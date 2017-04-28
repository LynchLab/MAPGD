#include "relatedness.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

///Takes a thing and does something.
size_t 
freqtoi(float_t in)
{
	return size_t(in*E_LIM*2) < E_LIM ? size_t(in*E_LIM*2) : E_LIM-1;
}

#ifdef EIGEN
void 
newton (Relatedness &a, std::map <Genotype_pair_tuple, size_t> &counts)
{
	Eigen::MatrixXd J(7,7);
	Eigen::MatrixXd iJ(7,7);
	Eigen::VectorXd R(7);

	std::map <Genotype_pair_tuple, size_t>::iterator it, end;
	Genotype_pair v;
	size_t c;

	float_t det, sumlnL;

	while (true){
		it=counts.begin();
		end=counts.end();
		J=Eigen::MatrixXd::Zero(7,7);
	while (it!=end ){
		std::cout << "[" << R << "]" << std::endl;
		v=Genotype_pair::from_tuple(it->first);
		c=it->second;
		J(0,0)+=J00(v, a)*c; J(0,1)+=J01(v, a)*c; J(0,2)+=J02(v, a)*c; J(0,3)+=J03(v, a)*c; 
				J(0,4)+=J04(v, a)*c; J(0,5)+=J05(v, a)*c; J(0,6)+=J06(v, a)*c; 
	
		J(1,0)+=J10(v, a)*c; J(1,1)+=J11(v, a)*c; J(1,2)+=J12(v, a)*c; J(1,3)+=J13(v, a)*c; 
				J(1,4)+=J14(v, a)*c; J(1,5)+=J15(v, a)*c; J(1,6)+=J16(v, a)*c; 

		J(2,0)+=J20(v, a)*c; J(2,1)+=J21(v, a)*c; J(2,2)+=J22(v, a)*c; J(2,3)+=J23(v, a)*c; 
				J(2,4)+=J24(v, a)*c; J(2,5)+=J25(v, a)*c; J(2,6)+=J26(v, a)*c; 

		J(3,0)+=J30(v, a)*c; J(3,1)+=J31(v, a)*c; J(3,2)+=J32(v, a)*c; J(3,3)+=J33(v, a)*c; 
				J(3,4)+=J34(v, a)*c; J(3,5)+=J35(v, a)*c; J(3,6)+=J36(v, a)*c; 

		J(4,0)+=J40(v, a)*c; J(4,1)+=J41(v, a)*c; J(4,2)+=J42(v, a)*c; J(4,3)+=J43(v, a)*c; 
				J(4,4)+=J44(v, a)*c; J(4,5)+=J45(v, a)*c; J(4,6)+=J46(v, a)*c; 

		J(5,0)+=J50(v, a)*c; J(5,1)+=J51(v, a)*c; J(5,2)+=J52(v, a)*c; J(5,3)+=J53(v, a)*c; 
				J(5,4)+=J54(v, a)*c; J(5,5)+=J55(v, a)*c; J(5,6)+=J56(v, a)*c; 

		J(6,0)+=J60(v, a)*c; J(6,1)+=J61(v, a)*c; J(6,2)+=J62(v, a)*c; J(6,3)+=J63(v, a)*c; 
				J(6,4)+=J64(v, a)*c; J(6,5)+=J65(v, a)*c; J(6,6)+=J66(v, a)*c; 

		R(0)+=H0(v, a);
		R(1)+=H1(v, a);
		R(2)+=H2(v, a);
		R(3)+=H3(v, a);
		R(4)+=H4(v, a);
		R(5)+=H5(v, a);
		R(6)+=H6(v, a);

		sumlnL+=lnL_NR(v, a);

		++it;
	};

	iJ=J.inverse();
	det=J.determinant();
	iJ/=det;
	R=iJ*R;

/*
	if (fabs(R[0])>B)
	{
		if(R[0]>0) R(0)=B;
			else R(0)=-B;
	}*/
	//params=[F_X, F_Y, T, g_XY, g_YX, d, D]

	a.f_X_-=R(0);
	a.f_Y_-=R(1);
	a.theta_XY_-=R(2);
	a.gamma_XY_-=R(3);
	a.gamma_YX_-=R(4);
	a.delta_XY_-=R(5);
	a.Delta_XY_-=R(6);
	}
}
#endif

std::map <Genotype_pair_tuple, size_t> 
hash_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y)
{
	std::stringstream fb_copy(file_buffer.str() );
	
	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	Population genotypes=gcf_in.read_header();
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		if (genotypes.likelihoods[x].N>1 && genotypes.likelihoods[y].N>1 ){
			counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 4)]+=1;
		}
	}
	return counts;
}

std::map <Genotype_pair_tuple, size_t> 
downsample_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y)
{
	std::stringstream fb_copy(file_buffer.str() );
	
	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	Population genotypes=gcf_in.read_header();
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		if (genotypes.likelihoods[x].N>1 && genotypes.likelihoods[y].N>1 ){
			counts[downvert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 2)]+=1;
		}
	}
	return counts;
}

Relatedness global_relatedness;

/*Does a regression of allele frequency of the samples on the population allele frequency*/
void 
set_e(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
	std::map<Genotype_pair_tuple, size_t>::iterator start=hashed_genotypes.begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=hashed_genotypes.end();
	std::map<Genotype_pair_tuple, size_t>::iterator it=start;

	Genotype_pair pair;

/*
        float_t X_MM;
        float_t X_Mm;
        float_t X_mm;
        float_t Y_MM;
        float_t Y_Mm;
        float_t Y_mm;
        float_t m;
*/
	float_t freq[E_LIM]={0}, sum[E_LIM]={0};
        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
		freq[freqtoi(pair.m)]+=pair.m;
		sum[freqtoi(pair.m)]+=it->second;
		relatedness.e_X_[freqtoi(pair.m)]+=exp(-pair.X_mm)+exp(-pair.X_Mm)/2.;
		relatedness.e_Y_[freqtoi(pair.m)]+=exp(-pair.Y_mm)+exp(-pair.Y_Mm)/2.;
		++it;
	}

	for (size_t x=0; x<E_LIM; ++x){
		relatedness.e_X_[x]=(relatedness.e_X_[x]-freq[x])/sum[x];
		relatedness.e_Y_[x]=(relatedness.e_Y_[x]-freq[x])/sum[x];
//		std::cerr << x << ": " << relatedness.e_X_[x] << ", " << relatedness.e_Y_[x] << std::endl;
	}
	global_relatedness=relatedness;
}


/*Guess starting values of relatedness for the maximization procedure*/
void 
gestimate(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &counts)
{
	relatedness.zero();
/*	std::map<Genotype_pair_tuple, size_t>::iterator start=counts.begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=counts.end();
	std::map<Genotype_pair_tuple, size_t>::iterator it=start;
	Genotype_pair pair;

	for (size_t x=0; x<?; x++){
	?
	double OX=
	double OY=

	double k1x=pow(OX*(1-OX), 2);
	double k2x=OX*(1-OX)*(3*OX*OX-3*OX+1);
	double k1y=pow(OY*(1-OY), 2);
	double k2y=OY*(1-OY)*(3*OY*OY-3*OY+1);

	k1=sqrt(k1x*k1y);
	k2=sqrt(k2x*k2y);

	gsl_matrix_set (Xm, c2, 0, k1);
	gsl_matrix_set (Xm, c2, 1, k2);
	gsl_vector_set (yv, c2, mu);
	gsl_vector_set (wv, c2, Den);

	}
#define beta(i) (gsl_vector_get(cv,(i)))
	gsl_multifit_wlinear (Xm, wv, yv, cv, cov, &chisq, work);
        std::cerr << "\t" << f_X_w << "\t" << f_Y/f_Y_w << "\t" << f_Y_w << "\t" << Theta/Theta_w << "\t" <<Theta_w <<"\t" << gamma_XY/gamma_XY_w << "\t" << gamma_XY_w << "\t" << gamma_YX/gamma_YX_w << "\t" << gamma_YX_w << "\t" << beta(1) << "\t" << "0" << "\t" << beta(0) << "\t" << std::endl;
*/
/*        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
               	inc_guess(relatedness, pair, it->second);
                ++(it);
                ++(it);
        }
*/

}

void 
inc_f(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	float_t P2=P*P;
	float_t var=P-P2;
	float_t denom=pow(var, 2);
	if(pair.m!=0){
		rel.f_X_+=(2.*(exp(-pair.X_Mm)/4.+exp(-pair.X_MM)-P2) -var)/denom*count;
		rel.f_Y_+=(2.*(exp(-pair.Y_Mm)/4.+exp(-pair.Y_MM)-P2) -var)/denom*count;
	};
}

void 
inc_theta(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	if(P!=0){
		rel.theta_XY_+=( (exp(-pair.X_Mm-pair.Y_Mm)/4.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,2) )/pow(P*(1-P),2)*count;
	};
}

void 
inc_gamma(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	if(P!=0 && P!=0.5){
		rel.gamma_XY_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/4.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1-P)*(rel.theta_XY_*(1+2*P)+rel.f_X_*P+P/(1-P) ) )/(P*(1-P) )/(0.5-P)/2*count;
		rel.gamma_YX_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/4.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1.-P)*(rel.theta_XY_*(1.+2.*P)+rel.f_Y_*P+P/(1.-P) ) )/(P*(1.-P) )/(0.5-P)/2.*count;
	};
}

void 
inc_Delta (Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{

}

//TODO FIX THIS TERRIBLE HACK JOB YOU SCHMUCK!

//int eval;

float_t
get_ll (const Relatedness &rel, const Genotype_pair &pair, const float_t count) 
{
//	eval++;
	/*This is the basic likelihood model that we are fitting. It essentially calculates the correlation coefficents
	for the first four moments of the joint distribution. The math needs to be cleaned up here. For now it is typed
	up to minimize the chance of typos. a, b, c and d are r.v. Representing the presence (1) or absence of (0) of
	the Major allele in the haploid genomes of individuals A (for a and b) and C (c and d). A is the r.v. defined as
	A~(a+b)/2 and C~(c+d)/2. Generally what we are doing here is calculating the expectations (e.g. E_A2) of the 
	joint distributions and using that to calculate the joint distribution itself. Problems can occur when 
	correlation coefficients are less than zero, because the probabilities of various observations (e.g. mm1mm2) can 
	become negative. These probabilities are forced to be zero, which may be a little arbitrary, but it seems to 
	work.*/


	float_t P, mm1mm2, Mm1mm2, MM1mm2, mm1Mm2, Mm1Mm2, MM1Mm2, mm1MM2, Mm1MM2, MM1MM2;

	P=1-pair.m;

	if (P==0) return 0;	
	float_t e=0;//global_relatedness.e_Y_[freqtoi(pair.m)]-global_relatedness.e_X_[freqtoi(pair.m)];
	//P+=(global_relatedness.e_Y_[freqtoi(pair.m)]+global_relatedness.e_X_[freqtoi(pair.m)]);

	/*This comes from the inverse matrix of the one used to calculate the moments.*/
	
	float_t A=P+e;//+e(P)*P; 	//mean major allele frequency in A
	float_t C=P-e;//-e(P)*P; 	//mean major allele frequency in C
//	std::cerr << freqtoi(pair.m) << ", " << e << ", " << P << ", " << A << ", " << C << std::endl;	
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
	std::vector < std::pair<Genotype_pair_tuple, size_t> > *hashed_genotypes_p=(std::vector <std::pair <Genotype_pair_tuple, size_t> > *) void_hashed_genotypes_p;
	Genotype_pair first;
	size_t count;

	rel.f_X_ = gsl_vector_get(v, 0);
	rel.f_Y_ = gsl_vector_get(v, 1);
	rel.theta_XY_ = gsl_vector_get(v, 2);
	rel.gamma_XY_ = gsl_vector_get(v, 3);
	rel.gamma_YX_ = gsl_vector_get(v, 4);
	rel.Delta_XY_ = gsl_vector_get(v, 5);
	rel.delta_XY_ = gsl_vector_get(v, 6);
	
	float_t sum=0;
	
	std::pair<Genotype_pair_tuple, size_t> *pair;
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator end=hashed_genotypes_p->end();
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator it=hashed_genotypes_p->begin();

//	#pragma omp parallel for private(first, count, pair) reduction(+:sum)
//	for (size_t x=0; x<hashed_genotypes_p->size(); x++){
	while(it!=end){
//		pair=&(*hashed_genotypes_p)[x];
//		pair=&it;
		first=Genotype_pair::from_tuple(it->first);
		count=it->second;
		sum+=get_ll(rel, first, count);
		it++;
	}
	if (std::isnan(sum) ) return FLT_MAX;
	return double(-sum);
}

/*Maximizes the relatedness*/
void 
maximize(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
	const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2rand;
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function gsl_func;
	
	size_t iter = 0;
	int status;
	double size;

	/* Starting point and stepsizes. I would really prefer to do this with a 
	 * Newton-Raphson method, which just inexplicably like more than the 
	 * Nelder-Mead, but I'm being lazy today, or more accurately there are 
	 * fundamental problems with setting the problem up to use the NR method.
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
	gsl_func.params=&hashed_genotypes_vector;

	s = gsl_multimin_fminimizer_alloc (T, 7);
	gsl_multimin_fminimizer_set (s, &gsl_func, x, ss);
	
//	eval=0;

	do {

		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
      
		if (status) break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-3);

	}  while (status == GSL_CONTINUE && iter < 800);

//	std::cerr << " Solved in " << eval << " evaluations " << std::endl;

	rel.f_X_ = gsl_vector_get(s->x, 0);
	rel.f_Y_ = gsl_vector_get(s->x, 1);
	rel.theta_XY_ = gsl_vector_get(s->x, 2);
	rel.gamma_XY_ = gsl_vector_get(s->x, 3);
	rel.gamma_YX_ = gsl_vector_get(s->x, 4);
	rel.Delta_XY_ = gsl_vector_get(s->x, 5);
	rel.delta_XY_ = gsl_vector_get(s->x, 6);
	
	rel.max_ll_=rel_ll(x, &hashed_genotypes_vector);
}

class small_rel{

	public:
	size_t skip;
	small_rel(size_t skip_) {
		skip=skip_;
	}
	double
	rel_ll_2 (const gsl_vector *v, void *void_hashed_genotypes_p)
	{
       		gsl_vector *x=gsl_vector_alloc(7);

		for (size_t y=0; y<7;y++){
			if(y==skip){
   		     		gsl_vector_set(x, y, 0);
			} else {
				gsl_vector_set(x, y, gsl_vector_get(v, y-(y<=skip) ) );
			}
		}
		return rel_ll(v, void_hashed_genotypes_p);
	}
};

void
get_llr(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> hashed_genotypes)
{
	const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function gsl_func;
	
	size_t iter = 0;
	int status;
	double size;

	std::vector <float_t> values={rel.f_X_, rel.f_Y_, rel.theta_XY_,rel.gamma_XY_,rel.gamma_YX_, rel.Delta_XY_, rel.delta_XY_};
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );


	/* Starting point and stepsizes. I would really prefer to do this with a 
	 * Newton-Raphson method, which just inexplicably like more than the 
	 * Nelder-Mead, but I'm being lazy today, or more accurately there are 
	 * fundamental problems with setting the problem up to use the NR method.
	 */
	for (size_t w=0; w<7; ++w) {

//		small_rel this_rel(w);
		
//		ss=gsl_vector_alloc(7);
		x=gsl_vector_alloc(7);

		rel.null_ll_ = rel_ll(x, &hashed_genotypes_vector);

	        gsl_vector_set(x, 0, rel.f_X_);
	        gsl_vector_set(x, 1, rel.f_Y_);
	        gsl_vector_set(x, 2, rel.theta_XY_);
	        gsl_vector_set(x, 3, rel.gamma_XY_);
	        gsl_vector_set(x, 4, rel.gamma_YX_);
	        gsl_vector_set(x, 5, rel.Delta_XY_);
	        gsl_vector_set(x, 6, rel.delta_XY_);

		rel.max_ll_ = rel_ll(x, &hashed_genotypes_vector);


//		gsl_vector_set_all(ss,0.15);	

//TODO FIX ME!
/*		gsl_func.n=6;
		gsl_func.f = rel_ll;
		gsl_func.params=&hashed_genotypes;

		s = gsl_multimin_fminimizer_alloc (T, 6);
		gsl_multimin_fminimizer_set (s, &gsl_func, x, ss);
		do {
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
      
			if (status) break;

			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, 1e-4);

		}  while (status == GSL_CONTINUE && iter < 600);
*/

		switch (w) {
			case 0:
	        		gsl_vector_set(x, 0, 0);
				rel.f_X_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.f_X_);
			break;
			case 1:
	        		gsl_vector_set(x, 1, 0);
				rel.f_Y_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.f_Y_);
			break;
			case 2:
	        		gsl_vector_set(x, 2, 0);
				rel.theta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.theta_XY_);
			break;
			case 3:
	        		gsl_vector_set(x, 3, 0);
				rel.gamma_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.gamma_XY_);
			break;
			case 4:
	        		gsl_vector_set(x, 4, 0);
				rel.gamma_YX_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.gamma_YX_);
			break;
			case 5:
	        		gsl_vector_set(x, 5, 0);
				rel.Delta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.Delta_XY_);
			break;
			case 6:
	        		gsl_vector_set(x, 6, 0);
				rel.delta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.delta_XY_);
			break;
		}
	}
	
}

#ifdef NOMPI
int estimateRel(int argc, char *argv[])
{
	std::cerr << "Warning: this program should generate AICc's. However, it doesn't and its _ll variables don't mean that much.\n"; 

	/* All the variables that can be set from the command line */

	std::string gcf_name="", rel_name="";

	Environment env;
	env.set_name("mapgd relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	gcf_name, "an error occurred while displaying the help message.", "input filename (default stdout)");
	env.optional_arg('o', "output", rel_name, "an error occurred while displaying the help message.", "output filename (default stdin)");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;		//The class to read to
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
	std::map <Genotype_pair_tuple, size_t> down_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();

	for (size_t x=0; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
			relatedness.set_X_name(x);
			relatedness.set_Y_name(y);
			hashed_genotypes=hash_genotypes(file_buffer, x, y);
			down_genotypes=downsample_genotypes(file_buffer, x, y);
			relatedness.zero();
			set_e(relatedness, hashed_genotypes);
		//	gestimate(relatedness, hashed_genotypes);
#ifdef EIGEN
			newton(relatedness, down_genotypes);
#else
			maximize(relatedness, down_genotypes);
#endif
		//	maximize(relatedness, hashed_genotypes);
			get_llr(relatedness, hashed_genotypes);

			rel_out.write(relatedness);
		}
	}
	rel_out.close();
	return 0;					//Since everything worked, return 0!.
}
#endif

#else 

int 
estimateRel(int argc, char *argv[])
{
	std::cerr << "This command depends on gsl, which could not be found. Please e-mail matthew.s.ackerman@gmail.com for help.\n";
	return 0;
}

#endif 
