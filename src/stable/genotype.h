#ifndef GENOTYPE_H_
#define GENOTYPE_H_
/* \breif genotypic likelihoods. Arbitrary ploidy
 */
#include "quartet.h"
#include "typedef.h"
#include <string>

class genotype {

private:
	std::string name_;
	uint8_t ploidy_;
	quartet_t *labeled_sequence_;
	float_t *probabilities_;
public:
	genotype (uint8_t ploidy){
		ploidy_=ploidy;
		probabilities_=new float_t ();
	};
};
#endif
