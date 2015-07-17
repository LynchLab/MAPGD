/* \breif genotypic likelihoods. Arbitrary ploidy
 */

#include "genotype.h"

void genotype::set_probabilities(const population &p, const quartet_t &q){
	float_t lnM, lnm;
	for (size_t x; x<probabilities_.size(); ++x){
		lnM=?;
		lnm=?;
		probabilities_[x]=q[p.major]*lnM+q[p.minor]*lnm+log(p.get_genotypic_frequency(x) )+lnbionomial_coef(x,ploidy) ;
	}
}	

float_t genotype::get_probabilities(const size_t &x) const {
	return probabilities_[x];
}	
