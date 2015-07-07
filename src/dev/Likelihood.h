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

//typedef double float_t;
/*class model {
	public 
};*/

class lnmultinomial {
private:
	std::vector <float_t> lnfact_vector;
	float_t *lnp_;
	count_t size_;
public:	
	lnmultinomial (float_t*, count_t);		//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (count_t);			//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (void);			//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t*);				//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t, float_t, float_t, float_t);	//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set(void (*)(allele_stat const &, float_t *), allele_stat const &);

	float_t & lnprob(const count_t *) const;	//returns the probabiltiy of the multinomial distribution . . .
	float_t & lnprob_approx(const count_t *) const;	//returns the probabiltiy of the multinomial distribution . . .
	float_t & lnfact(const count_t &) const;	//returns the log factorial of the count type numbers in the array
	float_t & lnmultinomcoef(const count_t *) const;//returns the log factorial of the count type numbers in the array
};


std::vector <std::pair <count_t, float_t> > sort (const float_t *, const count_t &);
std::vector <std::pair <count_t, count_t> > sort (const count_t *, const count_t &);
#endif
