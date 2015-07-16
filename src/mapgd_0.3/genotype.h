#ifndef GENOTYPE_H_
#define GENOTYPE_H_
/* \breif genotypic likelihoods. Arbitrary ploidy
 */
#include "quartet.h"
#include "typedef.h"
#include "binomial.h"
#include "stdint.h"
#include <string>

class genotype {

private:
	std::string name_;
	uint8_t ploidy_;
	float_t *probabilities_;
	binomial bin;
public:
	genotype (uint8_t ploidy){
		ploidy_=ploidy;
		bin=binomial(0.5);
		probabilities_=new float_t [bin.binom_coef(uint32_t(ploidy), uint32_t(2) )];
	};

	~genotype (){
		delete [] probabilities_;
	};
};
#endif
