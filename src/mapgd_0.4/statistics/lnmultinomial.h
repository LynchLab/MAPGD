//General likelihood functions

#ifndef LNMULTINOMIAL_H_
#define LNMULTINOMIAL_H_

#include <map>
#include <string>
#include <vector>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <limits>	// std::numeric_limits

#include "../typedef.h"
#include "../datatypes/allele.h"

/* \breif information that should be passed to vcf files.
 */

//* This is pathalogical. Working to kill it.*/
/* \breif returns the log... 
 */
class lnmultinomial {
private:
	std::vector <float_t> lnfact_vector;
	float_t *lnp_;
	count_t size_;
public:	
	lnmultinomial (float_t*, const count_t&);	//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (const count_t&);			//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	lnmultinomial (void);				//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	~lnmultinomial (void);				//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .

	void set (float_t*);				//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set (float_t, float_t, float_t, float_t);	//!< creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .
	void set(void (*)(const allele_t&, float_t *), const allele_t&);

	float_t lnprob(const count_t*) ;		//!< returns the probabiltiy of the multinomial distribution . . .
	float_t lnprob_approx(const count_t*);		//!< returns the probabiltiy of the multinomial distribution . . .
	float_t lnfact(const count_t&);			//!< returns the log factorial of the count type numbers in the array
	float_t lnmultinomcoef(const count_t*);		//!< returns the log factorial of the count type numbers in the array

	lnmultinomial& operator=(const lnmultinomial& rhs);  //!< a comment.
};

std::vector <std::pair <count_t, float_t> > sort (const float_t *, const count_t &);
std::vector <std::pair <count_t, count_t> > sort (const count_t *, const count_t &);
#endif
