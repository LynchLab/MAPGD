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
#include <iostream>

/* \breif information that should be passed to vcf files.
 */

//* This is pathalogical. Working to kill it.*/
class allele_stat { 
private:
public:
	bool pooled;		//Infered from pooled or labeled sequencing?

	char delim;		
	std::string id0;	//
	count_t id1;		//
	count_t excluded;	//

	allele_stat();		//

	float_t freq;		//frequency of major allele.
	count_t minor;		//idenity of minor allele.
	count_t major;		//idenity major allele.
	float_t error;		//ml error rate.
	float_t null_error;	//ml error rate.
	float_t coverage;	//population coverage.
	float_t ll, monoll, hwell; //loglikelihoods.

	//FOR LABELED SEQUENCING ONLY!
	float_t MM; 		//frequency of genotype MM in the population.
	float_t Mm; 		//frequency of genotype Mm in the population.
	float_t mm; 		//frequency of genotype mm in the popualtion.
	float_t N;		//number of individual at the site.
	float_t f;		//HW statistic.
	float_t h;		//heterozygosity.
	float_t gof; 		//gof statistic.
	float_t efc; 		//number of 'effective' chromosomes in the sample.

	allele_stat & operator=(const allele_stat &);
	friend std::ostream& operator<< (std::ostream&, const allele_stat&);
};

/* \breif returns the log... 
 */
class lnmultinomial {
private:
	std::vector <float_t> lnfact_vector;
	float_t *lnp_;
	count_t size_;
public:	
	lnmultinomial (float_t*, const count_t&);	//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (const count_t&);			//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (void);				//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	~lnmultinomial (void);				//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .

	void set (float_t*);				//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t, float_t, float_t, float_t);	//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set(void (*)(const allele_stat&, float_t *), const allele_stat&);

	float_t lnprob(const count_t*) ;		//returns the probabiltiy of the multinomial distribution . . .
	float_t lnprob_approx(const count_t*);		//returns the probabiltiy of the multinomial distribution . . .
	float_t lnfact(const count_t&);			//returns the log factorial of the count type numbers in the array
	float_t lnmultinomcoef(const count_t*);		//returns the log factorial of the count type numbers in the array
};

std::vector <std::pair <count_t, float_t> > sort (const float_t *, const count_t &);
std::vector <std::pair <count_t, count_t> > sort (const count_t *, const count_t &);
#endif
