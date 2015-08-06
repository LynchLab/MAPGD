#include "allele.h"

allele_t & allele_t::operator=(const allele_t & rhs) {
	if (this != &rhs) { 
		memcpy(this, &rhs, sizeof(allele_t) );
	}
	return *this;
};
