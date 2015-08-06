#ifndef BINOMIAL_H_
#define BINOMIAL_H_

#include <vector>
#include "typedef.h"

/// A class that returns the binomial coefficent or the pdf the binomial distribution. Not Log.
/*!	\breif A class to return the binomial coefficent and binomial probabilities. NOT LOG!
 */
class binomial {
private:
	std::vector <uint32_t> fact_vector;			//!< The vector that stores look up values for the binomial coefficent.
	real_t p_;
	uint32_t *bp_;		
	uint32_t size_;						
public:	
	binomial () {};
	binomial (const real_t &p) {p_=p;};			//!< The constructor can be called with a real_t to specify the probability of success.
	~binomial (void){ fact_vector.clear(); }		

	uint32_t fact(const uint32_t&) ;			//!< Returns the probabiltiy of the binomial distribution . . .
	uint32_t binom_coef(const uint32_t&, const uint32_t&);	//!< Returns the binomial coefficent.
};

#endif
