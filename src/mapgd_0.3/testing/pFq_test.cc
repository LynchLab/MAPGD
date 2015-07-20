#include "../individual-likelihood.h"
#include <list>

//#include "gsl/gsl_math.h"		The GSL implementation of the hypergeometric function is 
//#include "gsl/gsl_sf_hyperg.h"	officially completely broken. It returns exceptions over nearly 
//					the entire parameter space of the function. I may end up using
//					the series expansion.... *sigh*


int main (int argc, char *argv[]){
	
	count_t l[4];
	float_t tP;
	float_t E=0, V=0;
	allele_stat a;
	models model;

	a.error=0.01;

	a.MM=atof(argv[1]);
	a.Mm=atof(argv[2]);
	a.mm=atof(argv[3]);

	count_t N_=atoi(argv[4]);

	a.major=a.MM+a.Mm/2.;
	a.minor=1-a.major;

	std::cout << pFq::_2F1(-2., -2., 1, 0.5) << std::endl; 
	std::cout << pFq::_2F1_simp(-2., 0.5) << std::endl; 
	
	float_t H=(0.5*(1.-a.error)+0.5*a.error/3.);
	std::cout <<  "E " << a.error << ", " << N_ << "\n" <<
		pow(a.MM,2)*pow(1.-a.error-1., 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(1-a.error, 2) /pow( (1.-a.error)-1.,2))
		+pow(a.mm,2)*pow(a.error-1, 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(a.error, 2) /pow(a.error-1.,2))
		+pow(a.Mm,2)*pow(H-1., 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(H, 2) /pow(H-1.,2))
		 << '\n';

	for (size_t x=0; x<N_+1; ++x){
		memset(l, 0, sizeof(count_t)*4);
		l[a.major]=x;
		l[a.minor]=N_-x;
		tP=model.lnP(l, a);
		std::cout << "x : " << x << ", " << exp(tP) << std::endl;
		E+=exp(tP)*exp(tP);
		V+=pow(exp(tP), 3);
	}
	V-=pow(E, 2);
	std::cout << "O " << E << ", " << V << "\n\n";
}

