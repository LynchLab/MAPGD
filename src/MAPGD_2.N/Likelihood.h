//General likelihood functions

#ifndef LIKELIHOOD_H_
#define LIKELIHOOD_H_

#include <map>
#include "math.h"
#include "typedef.h"
#include <string>
#include <vector>
#include <cfloat>
#include <algorithm>
//typedef double float_t;

/*class model {
	public 
};*/

typedef struct allele_stat { 
	bool pooled;		//Infered from pooled or labeled sequencing?

	float_t freq;		//frequency of major allele.
	count_t minor;		//idenity of minor allele.
	count_t major;		//idenity major allele.
	float_t error;		//ml error rate.
	float_t null_error;	//ml error rate.
	float_t coverage;	//population coverage.
	float_t ll; 		//loglikelihood.

	//FOR LABELED SEQUENCING ONLY!
	float_t MM; 		//frequency of genotype MM in the population.
	float_t Mm; 		//frequency of genotype Mm in the population.
	float_t mm; 		//frequency of genotype mm in the popualtion.
	float_t N; 		//number of individual at the site.
	float_t f; 		//HW statistic.
	float_t gof; 		//gof statistic.
	float_t efc; 		//number of 'effective' chromosomes in the sample.

} allele_stat_t;

class lnmultinomial {
private:
	std::vector <float_t> lnfact_vector;
	float_t *lnp;
	count_t _size;
public:	
	lnmultinomial (float_t*, count_t);		//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (count_t);			//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t*);				//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t, float_t, float_t, float_t);	//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set(void (*)(allele_stat const &, float_t *), allele_stat const &);
	float_t lnprob(const count_t *);		//returns the probabiltiy of the multinomial distribution . . .
	float_t lnprob_approx(const count_t *);		//returns the probabiltiy of the multinomial distribution . . .
	float_t lnfact(const count_t &);		//returns the log factorial of the count type numbers in the array
	float_t lnmultinomcoef(const count_t *);	//returns the log factorial of the count type numbers in the array
};

std::vector <std::pair <count_t, float_t> > sort (const float_t *, const count_t &);
std::vector <std::pair <count_t, count_t> > sort (const count_t *, const count_t &);
#endif
