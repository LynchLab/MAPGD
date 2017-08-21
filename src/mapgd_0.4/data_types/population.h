#ifndef _POPULATION_H_
#define _POPULATION_H_

#include <string.h>
#include <iostream>
#include <sstream>
#include "data.h"
#include "sample_name.h"
#include "typedef.h"
#include "genotype.h"

/// Population genotypes.
class Population : public Indexed_data {
private:
	std::vector <Genotype>::iterator igl_;		//!< an iterator to allow us to iterate over the likelihoods.
	std::vector <std::string> sample_names_;	//!< a vector of sample names.
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Population(Columns);
	}
public:
//Data block
	std::vector <Genotype> likelihoods;		//!< Genotypic likelihood

	gt_t major;					//!< identity of the major allele
	gt_t minor;					//!< identity of the minor allele
	float_t m;					//!< minor allele frequency
	float_t f;					//!< departure from HWE
//end Data block
	Population();					//!< simple constructor.
	Population(const Sample_name &);		//!< constructor using sample names. 
	Population(const std::vector <std::string> &);	//!< constructor needed by map_file. String should be coloumn names. 
	Population(const Population &); 		//!< constructor using a Population_Genotype
	~Population();					//!< destructor.
	size_t size() const;					//!< Returns the number of samples.
	void add(const Genotype &likelihood);			//!< append a sample to the likelihood.
	void add(const float_t &lMM, const float_t &lMm, const float_t &lmm, const count_t &N); //!< append a sample to the likelihood.
	void clear();						//!< clear likelihoods.
	std::string header(void) const;				//!< print header.

	static const std::string table_name;			//!< destination table in Db.
	static const std::string file_name;			//!< defualt file extention.
	static const bool binary;
	const bool get_binary() const;

	inline std::vector <std::string> get_sample_names(void) const {return sample_names_;};		//!< names of the samples sequenced.
	inline void set_sample_names(const std::vector <std::string>& sample_names) {
		sample_names_=sample_names;
		likelihoods.resize(sample_names.size() );
	};		//!< names of the samples sequenced.

	Population & operator= (const Population&);

	void write (std::ostream&) const;
	void read (std::istream&);

	const std::string get_file_name() const {return Population::file_name;};
        const std::string get_table_name() const {return Population::table_name;};

	void write_binary (std::ostream&) const;
	void read_binary (std::istream&);
};

#endif
