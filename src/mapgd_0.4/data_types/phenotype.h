#ifndef _PHENOTYPE_H_
#define _PHENOTYPE_H_

#include <cstring>

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#ifdef EIGEN
//#include "Eigen/Core"
#endif

#include "typedef.h"
#include "data.h"
#include "stream_tools.h"

/// Phenotype data.
class Phenotype : public Data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Phenotype(Columns);
	}
	size_t n_samples_, n_traits_;
public:
	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

	std::vector <std::string> sample_name;
	std::vector <std::string> trait;

#ifdef EIGEN
	//Eigen::MatrixXf value;
#else
	std::vector <std::vector <real_t > > value;
#endif

	Phenotype();	
	Phenotype(const std::vector <std::string> &); 
	Phenotype(const size_t &); 

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;


	void add_sample(const uint32_t &, const real_t *);
	//! zeros values and sets names to empty.
//	void clear(void); 
	//! zeros values, but doesn't set names to empty.
//	void zero(void);  

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;

	static const bool binary;

	const bool get_binary() const;

//	const size_t sample_size(void) const;

};

#endif
