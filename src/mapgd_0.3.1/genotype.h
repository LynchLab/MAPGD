#ifndef GENOTYPE_H_
#define GENOTYPE_H_
#include "quartet.h"
#include "typedef.h"
#include "binomial.h"
#include "stdint.h"
#include "allele_stat.h"
#include <string>


/*!	\breif A data type to store genotypic probabilities. 
 *
 *	For members to be stored in map files they must contain only fixed width . . . Individuals can have names that appear in the column header. 
 */
class genotype {

private:
	std::string name_;
	uint8_t ploidy_;
	std::vector<float_t> probabilities_;
	//binomial bin;
public:
	/*!	\breif The constructor requries the ploidy of the individual.
	 */
	genotype (uint8_t ploidy)
	{
		ploidy_=ploidy;
		probabilities_=std::vector <float_t>(ploidy_+1);
	};

	~genotype () { };

	void set_probabilities(const allele_stat &, const quartet_t &);		//!< Set genotypic probabilities. Needs a likelihood function, prior genotypic probabilites, and a quartet. 
	float_t get_probabilities(const size_t &x) const;			//!< Returns genotypic probabilities. 
};
#endif
