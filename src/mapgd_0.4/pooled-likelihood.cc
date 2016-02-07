#include "pooled-likelihood.h"

	
/*! \breif sets the probabilities for reads in pooled sequence at a polymorphic site.
*/

void polymorphicmodel(Allele const &a, float_t *l)
{
	float_t e3=(a.error/3.);
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=( (1.-a.error)*a.freq+e3*(1.-a.freq) );
	l[a.minor]=( (1.-a.freq)*(1.-a.error)+e3*a.freq );
	l[0]=log(l[0]);
	l[1]=log(l[1]);
	l[2]=log(l[2]);
	l[3]=log(l[3]);
}

/*! \breif sets the probabilities for reads in pooled sequence at a site fixed for the major allele.
*/
void monomorphicmodel(Allele const &a, float_t *l)
{
	float_t e3=(a.error/3.);
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=(1.-a.error);
	l[0]=log(l[0]);
	l[1]=log(l[1]);
	l[2]=log(l[2]);
	l[3]=log(l[3]);
}

/*! \breif sets the probabilities for reads in pooled sequence at a site fixed for the minor allele.
*/
void fixedmorphicmodel(Allele const &a, float_t *l)
{
	float_t e3=(a.error/3.);
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.minor]=(1.-a.error);
	l[0]=log(l[0]);
	l[1]=log(l[1]);
	l[2]=log(l[2]);
	l[3]=log(l[3]);
}

