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
	genotype & operator= (const genotype&);
};

std::ostream& operator<< (std::ostream& out, const genotype& x);
std::istream& operator>> (std::istream& in, genotype& x);

class population_genotypes{
private:
	std::vector <genotype>::iterator igl_;		//!< an iterator to allow us to iterate over the likelihoods.
	std::vector <std::string> sample_names_;	//!< a vector of sample names.
	bool frozen;					//!< a flag to indicate that no more samples will be added 
public:
	std::vector <genotype> likelihoods;		//!< genotypic likelihood
	std::vector <std::string> sample_names;		//!< genotypic likelihood

	gt_t major;					//!< identity of the major allele
	gt_t minor;					//!< identity of the minor allele
	float_t m;					//!< minor allele frequency
	id0_t id0;					//!< scaffold number
	id1_t id1;					//!< location on scaffold

	population_genotypes();					//!< simple constructor.
	population_genotypes(std::vector <std::string>);	//!< constructor needed by map-file. String should be coloumn names. 
	population_genotypes(const population_genotypes &); 	//!< constructor using a population_genotype
	~population_genotypes();				//!< destructor.
	size_t size() const;					//!< Returns the number of samples.
	void add(const genotype &likelihood);			//!< append a sample to the likelihood.
	void add(const float_t &lMM, const float_t &lMm, const float_t &lmm, const count_t &N); //!< append a sample to the likelihood.
	void clear();						//!< clear likelihoods.
	std::string header(void) const;				//!< print header.
	static const std::string table_name;			//!< destination table in Db.
	static const std::string file_name;			//!< defualt file extention.
	inline std::vector <std::string> get_sample_names(void) const {return sample_names_;};		//!< names of the samples sequenced.
	inline void set_sample_names(const std::vector <std::string>& sample_names) {
		sample_names_=sample_names;
		likelihoods.resize(sample_names.size() );
	};		//!< names of the samples sequenced.
};

std::ostream& operator<< (std::ostream& out, const population_genotypes& x);
std::istream& operator>> (std::istream& in, population_genotypes& x);

/* Not using an enum to avoid a static cast later on*/
typedef std::tuple <float_t, float_t, float_t, float_t, float_t, float_t, float_t> Genotype_pair_tuple; 

class Genotype_pair {
public:
	Genotype_pair(){};
	Genotype_pair(const float_t&, const float_t &, const float_t &, const float_t &, const float_t &, const float_t &, const float_t &);
	float_t X_MM;
	float_t X_Mm;
	float_t X_mm;
	float_t Y_MM;
	float_t Y_Mm;
	float_t Y_mm;
	float_t m;
	static Genotype_pair_tuple to_tuple(const Genotype_pair &);
	static Genotype_pair from_tuple(const Genotype_pair_tuple &);
};


Genotype_pair_tuple convert(const genotype &, const genotype &, const float_t &, const uint8_t &);
#endif  
