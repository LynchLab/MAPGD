/* synonym for population? */

#ifndef SAMPLE_GOF_H
#define SAMPLE_GOF_H

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "data.h"

///	The sample specific goodness of fit.
class Sample_gof : public Data { 
private:
	void write (std::ostream&) const;	//!< use the << operator to write Allele.
	void read (std::istream&);		//!< use the >> operator to read Allele.
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Sample_gof(Columns);
	}
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	float_t number_;
	std::string name_;

	Sample_gof();	
	Sample_gof(std::vector <std::string>) : Sample_gof(){};	
	Sample_gof(const std::string &, const float_t &);	

	std::string header(void) const;
	size_t size(void) const;

	static const std::string file_name;					//!< The dafualt extention for files.
	static const std::string table_name;					//!< Destination table in Db.
	static const bool binary;					//!< Destination table in Db.

	const bool get_binary () const;					//!< Destination table in Db.
};

#endif
