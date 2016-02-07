/* synonym for population? */

#ifndef _POOLED_DATA_H_
#define _POOLED_DATA_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "base.h"
#include "data_types/allele.h"

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class Pooled_data { 
private:
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	Base major, minor; 
	id0_t id0;
	id1_t id1;
	count_t coverage;
	float_t error;

	std::vector <std::string> names_;
	std::vector <float_t> p;
	std::vector <float_t> polyll;
	std::vector <float_t> majorll;
	std::vector <float_t> minorll;
	
	Pooled_data();	
	Pooled_data(const std::vector <std::string> &);
	
	Allele to_allele(const size_t &);

	void set_sample_names(const std::vector <std::string> &);

	std::string header(void) const;
	size_t size(void) const;

	friend std::ostream& operator<< (std::ostream&, const Pooled_data&);	//!< use the << operator to write Allele.
	friend std::istream& operator>> (std::istream&, Pooled_data&);		//!< use the >> operator to read Allele.
	static const std::string file_name;					//!< The dafualt extention for files.
	static const std::string table_name;					//!< Destination table in Db.
};

#endif
