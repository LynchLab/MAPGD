#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_

#include "typedef.h"
#include <iostream>

class Genotype {
private :
public :
	float_t MM, Mm, mm;						//!< Major Major, Major minor, minor minor
	count_t N;							//!< total depth of coverage.
	Genotype(const float_t &MM, const float_t &Mm, const float_t &mm, const count_t &N);	//!< constructor.
	Genotype();							//!< constructor.
	Genotype & operator= (const Genotype&);
	friend std::ostream& operator << (std::ostream& out, const Genotype& x);
	friend std::istream& operator >> (std::istream& in, Genotype& x);
};

#endif
