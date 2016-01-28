#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_ 

#include <list>
#include <vector>
#include <tuple>
#include "typedef.h"

/** @brief genotypic likelihoods.
  *
  */

class genotype{
private :
public :
	float_t lMM, lMm, lmm;						//!< Major Major, Major minor, minor minor
	count_t N;							//!< total depth of coverage.
	genotype(const float_t &MM, const float_t &Mm, const float_t &mm, const count_t &N);	//!< constructor.
	genotype();							//!< constructor.
};
std::ostream& operator<< (std::ostream& out, const genotype& x);
std::istream& operator>> (std::istream& in, genotype& x);

class population_genotypes{
private:
	std::vector <genotype>::iterator igl;	//!<
	std::vector <std::string> sample_names_;
	bool frozen;				//!< ?
public:
	std::vector <genotype> likelihoods;		//!< genotypic likelihood
	std::vector <std::string> sample_names;		//!< genotypic likelihood

	count_t major;
	count_t minor;
	float_t m;					//!< minor allele frequency
	id0_t id0;
	id1_t id1;

	population_genotypes();				//!< constructor.
	population_genotypes(std::vector <std::string>);
	population_genotypes(const population_genotypes & popgl);		//!< constructor.
	~population_genotypes();				//!< destructor.
	size_t size() const;				//!< 
	void add(const genotype &likelihood);			//!< append a likelihood.
	void add(const float_t &lMM, const float_t &lMm, const float_t &lmm, const count_t &N);
	void clear();
	std::string header(void) const;
	static const std::string table_name;
	static const std::string file_name;
	inline std::vector <std::string> get_sample_names(void) const {return sample_names_;};		//!< names of the samples sequenced.
	inline void set_sample_names(const std::vector <std::string>& sample_names) {
		sample_names_=sample_names;
		likelihoods.resize(sample_names.size() );
	};		//!< names of the samples sequenced.
};

std::ostream& operator<< (std::ostream& out, const population_genotypes& x);
std::istream& operator>> (std::istream& in, population_genotypes& x);



typedef std::tuple <float_t, float_t, float_t, float_t, float_t, float_t, float_t> genotype_pair; 
genotype_pair convert(population_genotypes&, const count_t &, const count_t &, const float_t &, const float_t &, const float_t &, const float_t &ll_t);
#endif  
