#include "pooledLikelihood.h"
		
void polymorphicmodel(allele_stat const &a, float_t *l){
	//TODO FIX THIS!
	float_t e3=a.error/3.;
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=(1.-a.error)*a.freq+a.error*(1.-a.freq);
	l[a.minor]=(1.-a.freq)*(1.-a.error)+a.error*a.freq;
};

void monomorphicmodel(allele_stat const &a, float_t *l){
	//TODO FIX THIS!
	float_t e3=a.error/3.;
	l[0]=e3;
	l[1]=e3;
	l[2]=e3;
	l[3]=e3;
	l[a.major]=(1.-a.error);
};

