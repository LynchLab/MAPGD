#ifndef GENOTYPE_H_
#define GENOTYPE_H_

#include <string>

#include "quartet.h"
#include "allele.h"

#include "../typedef.h"
#include "../statistics/binomial.h"

/*!	\breif A data type to store genotypic probabilities. 
 *
 *	For members to be stored in map files they must contain only fixed width . . . Individuals can have names that appear in the column header. 
 */
class genotype {

private:
	std::string name_;
	uint8_t ploidy_;
	std::vector<real_t> probabilities_;
	//binomial bin;
public:
	/*!	\breif The constructor requries the ploidy of the individual.
	 */
	genotype (uint8_t ploidy)
	{
		ploidy_=ploidy;
		probabilities_=std::vector <real_t>(ploidy_+1);
	};

	~genotype () { };

	void set_probabilities(const allele &, const quartet_t &);		//!< Set genotypic probabilities. Needs a likelihood function, prior genotypic probabilites, and a quartet. 
	real_t get_probabilities(const size_t &x) const;			//!< Returns genotypic probabilities. 
};
#endif
