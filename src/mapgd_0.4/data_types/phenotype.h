#ifndef _PHENOTYPE_H_
#define _PHENOTYPE_H_

#include <string.h>
#include <iostream>
#include <sstream>
#include "data.h"
#include "typedef.h"

/// Population phenotypes.
class Phenotype : public Data {
private:
	std::vector <std::string> sample_names_;	//!< a vector of sample names.
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Population(Columns);
	}
public:
	std::vector <real_t> z_raw;	//!< uncentered phenotypes
	std::vector <real_t> z;	//!< centered phenotypes
	std::vector <real_t> z_prime;	//!< centered w/ inbreeding removed
	std::vector <real_t> a_hat;	//!< additive genetic value
	std::vector <real_t> d_hat;	//!< dominance deviation
	std::vector <real_t> e_hat;	//!< environmental deviation

	Phenotype();					//!< simple constructor.
	Phenotype(const std::vector <std::string> &);	//!< constructor needed by map-file. String should be coloumn names. 
	Phenotype(const Phenotype &); 		//!< constructor using a Phenotypes
	~Population();					//!< destructor.
	size_t size() const;					//!< Returns the number of samples.
	std::string header(void) const;				//!< print header.

	static const std::string table_name;			//!< destination table in Db.
	static const std::string file_name;			//!< defualt file extention.
	static const bool binary;

	const bool get_binary() const;

	inline std::vector <std::string> get_sample_names(void) const {return sample_names_;};		//!< names of the samples sequenced.
	inline void set_sample_names(const std::vector <std::string>& sample_names) {
		sample_names_=sample_names;
		likelihoods.resize(sample_names.size() );
	};		

	Phenotype & operator= (const Phenotype&);
	void write (std::ostream&) const;
	void read (std::istream&);
};

#endif
