#include "pooled_likelihood.h"

	
/*! \breif sets the probabilities for reads in pooled sequence at a polymorphic site.
*/

void polymorphicmodel(allele_t const &a, float_t *l){
	float_t e3=a.error/3.;
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=(1.-a.error)*get_freq(a)+e3*(1.-get_freq(a) );
	l[a.minor]=(1.-get_freq(a) )*(1.-a.error)+e3*get_freq(a);
};

/*! \breif sets the probabilities for reads in pooled sequence at a site fixed for the major allele_t.
*/
void monomorphicmodel(allele_t const &a, float_t *l){
	float_t e3=a.error/3.;
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=(1.-a.error);
};

/*! \breif sets the probabilities for reads in pooled sequence at a site fixed for the minor allele_t.
*/
void fixedmorphicmodel(allele_t const &a, float_t *l){
	float_t e3=a.error/3.;
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.minor]=(1.-a.error);
};

