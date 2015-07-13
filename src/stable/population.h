#ifndef POPULATION_H_
#define POPULATION_H_

#include "genotype.h"
#include <vector>
#include <string>

class population {

private:
	std::string name_;
	uint16_t population_size_;
	std::vector <const genotype*> genotypes_;
	const quartet_t* pooled_sequence_;
public:
	population (const quartet_t&);
	population (const count_t&);
};
#endif
