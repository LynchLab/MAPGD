#ifndef MODELS_H_
#define MODELS_H_
/*! \breif this class perfoms some likelihood calculations. It keeps some look-up tables that will grow as calculations
 *	are performed, and may become to large over time. TODO: implement automatic look-up table purging. 
 *
 *  Currently likelihood calculations are all performed by a set of 'models' that we give to the "multinomial" class via the 
 *  multinomial::set method. These models should be able to look at the allele_stat structure, which contains information about
 *  the error rate at the locus and the identity of the major and minor allele, and return a set of four log probabilities of 
 *  observing each particular nucleotide in a given call.
*/	

#include "lnmultinomial.h"
#include "allele_stat.h"
#include "typedef.h"
#include "quartet.h"
#include "locus.h"

class models{
private:
lnmultinomial lnMM_, lnMm_, lnmm_;		//!< The three multinomials we will use for probability calculations.
						// The '4' specifies the number of categories of the distribution.
						// Since these represent the distribution of the four nucleotides 
						// A, C, G and T, we use 4 categories.
lnmultinomial lnMMP_, lnMmP_, lnmmP_;		

float_t E0_, E1_, E2_;			//!< Values used in calculations;
public:

	models(void);
	~models(void);

	/*! \breif Gets the log likelihood of the observations at Locus, given allele_stat.*/
	float_t loglikelihood(const Locus &, const allele_stat &, const count_t &);

	/*! \breif Gets the log likelihood of the observations at Locus, given allele_stat.*/
	float_t goodness_of_fit (Locus &, const allele_stat &, std::vector <float_t> &, const count_t &, const float_t &);

	/*! \breif Initilizes the ? for goodness of fit calculations.*/
	void init_gof(const count_t *, const allele_stat &);

	/*! \breif Returns the likelihood of a ?? goodness of fit... blah blah blah.*/
	float_t get_gof(const count_t *, const allele_stat &);

	/*! \breif Clalculates genotypic likelihoods. Not implement, may be depricated.*/
	float_t genotypelikelihood(const quartet_t &, const allele_stat &, const count_t &);

	/*! \breif Clalculates goodness of fit likelihoods.*/
	inline float_t lnL(const float_t &logMM, const float_t &logMm, const float &logmm, const count_t *count){
		/*posterior = prior x likelihood */
		E0_=logMM+lnMMP_.lnprob(count);
		E1_=logMm+lnMmP_.lnprob(count);
		E2_=logmm+lnmmP_.lnprob(count);

		if (E0_>E2_) std::swap(E2_, E0_);
		if (E1_>E2_) std::swap(E2_, E1_);
					
		//We make E2 the largest (i.e. least negative) value of the three. E0 and E1 will often times be 
		//extreamly small, and can even be negative infinity. So we are basically just ensuring that we 
		//don't throw away a lot of precission by returning log(exp(E2) ). 

		return log(exp(E0_-E2_)+exp(E1_-E2_)+1.)+E2_;
	}

	models& operator=(const models& rhs);  
};

/*! \breif A functor to set probabilities in the lnmultinomial members of model: homozygous minor.*/
float_t *mmModel(allele_stat *);
/*! \breif A functor to set probabilities in the lnmultinomial members of model: heterozygous.*/
float_t *MmModel(allele_stat *);
/*! \breif A functor to set probabilities in the lnmultinomial members of model: homozygous major.*/
float_t *MMModel(allele_stat *);
/*! \breif A functor to set probabilities in the lnmultinomial members of model: homozygous minor.*/
float_t *mmModelP(allele_stat *);
/*! \breif A functor to set probabilities in the lnmultinomial members of model: heterozygous.*/
float_t *MmModelP(allele_stat *);
/*! \breif A functor to set probabilities in the lnmultinomial members of model: homozygous major.*/
float_t *MMModelP(allele_stat *);
#endif
