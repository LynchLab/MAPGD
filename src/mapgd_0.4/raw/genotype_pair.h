#ifndef _GENOTYPE_PAIR_H_
#define _GENOTYPE_PAIR_H_ 

#include <list>
#include <vector>
#include <tuple>
#include <sstream>

#include "typedef.h"
#include "genotype.h"

/** @brief genotypic likelihoods.
  *
  */

/* Not using an enum to avoid a static cast later on*/
typedef std::tuple <float_t, float_t, float_t, float_t, float_t, float_t, float_t> Genotype_pair_tuple; 

class Genotype_pair {
public:
	Genotype_pair(){};
	Genotype_pair(const float_t&, const float_t &, const float_t &, const float_t &, const float_t &, const float_t &, const float_t &);
	float_t X_MM;
	float_t X_Mm;
	float_t X_mm;
	float_t Y_MM;
	float_t Y_Mm;
	float_t Y_mm;
	float_t m;
	static Genotype_pair_tuple to_tuple(const Genotype_pair &);
	static Genotype_pair from_tuple(const Genotype_pair_tuple &);
};


Genotype_pair_tuple convert(const Genotype &, const Genotype &, const float_t &, const uint8_t &);
Genotype_pair_tuple downvert(const Genotype &, const Genotype &, const float_t &, const uint8_t &);
#endif  
