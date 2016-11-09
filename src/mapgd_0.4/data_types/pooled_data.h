#ifndef _POOLED_DATA_H_
#define _POOLED_DATA_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <math.h>

#include "typedef.h"
#include "base.h"
#include "data.h"
#include "allele.h"

///	Information genereated from pooled data.
/** This is information generated from pooled population data (i.e. mash everybody up together and sequence the slurry).
 */
class Pooled_data : public Indexed_data { 
private:
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Pooled_data(Columns);
	}
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	Base major, minor; 
	count_t coverage;
	float_t error;

	std::vector <std::string> names_;//!< the sample names.
	std::vector <float_t> p;	//!< a vector of allele frequencies
	std::vector <float_t> cov;	//!< a vector of coverage
	std::vector <float_t> polyll;	//!< a vector of polymorphic log likelihoods
	std::vector <float_t> fixedll;	//!< a vector of fixed log likelihoods 
	
	Pooled_data();	
	Pooled_data(const std::vector <std::string> &);
	
	Allele to_allele(const size_t &);

	void set_sample_names(const std::vector <std::string> &);
	Pooled_data & operator=(const Pooled_data &);	//!< use the = operator to assign Pooled_data.

	std::string header(void) const;
	size_t size(void) const;

	//! used to write Allele. Inherits << from Data
	void write (std::ostream&) const;	
	//! used to read Allele. Inherits >> from Data
	void read (std::istream&);
	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.
	static const bool binary;
	const bool get_binary() const;
};

#endif
