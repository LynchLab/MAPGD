//General likelihood functions

#ifndef _LNMULTINOMIAL_H_
#define _LNMULTINOMIAL_H_

#include <map>
#include <string>
#include <vector>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <string.h>
#include <limits>	// std::numeric_limits

#include "typedef.h"
#include "data_types/allele.h"

/* \brief information that should be passed to vcf files.
 */

//* This is pathological. Working to kill it.*/
/* \brief returns the log... 
 */
class lnmultinomial {
private:
	std::vector <float_t> lnfact_vector;		//!< A look up table for log factorial values.
	float_t *lnp_;					//!< TODO Add a description. 
	count_t size_;					//!< The number of categories in the multinomial.
public:	
	lnmultinomial (float_t*, const count_t&);	//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	lnmultinomial (const count_t&);			//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	lnmultinomial (void);				//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	~lnmultinomial (void);				//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	void set (float_t*);				//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	void set (float_t, float_t, float_t, float_t);	//!< Creates a function that returns log probabilities from a multinomial distribution with parameters float_t . . .
	void set(void (*)(const Allele&, float_t *), const Allele&);

	float_t lnprob(const count_t*) ;		//!< Returns the probability of the multinomial distribution.
	float_t lnprob_approx(const count_t*);		//!< Returns the probability of the multinomial distribution ().
	float_t lnfact(const count_t&);			//!< Returns the log factorial of the count type numbers in the array
	float_t lnmultinomcoef(const count_t*);		//!< Returns the log factorial of the count type numbers in the array

	lnmultinomial& operator=(const lnmultinomial& rhs);  //!< Use = operator to copy a lnmultinomial.
};

std::vector <std::pair <count_t, float_t> > sort (const float_t *, const count_t &);
std::vector <std::pair <count_t, count_t> > sort (const count_t *, const count_t &);
#endif
