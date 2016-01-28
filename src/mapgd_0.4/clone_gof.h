/* synonym for population? */

#ifndef CLONE_GOF_H
#define CLONE_GOF_H

#include <iostream>
#include "typedef.h"
#include "pro-file.h"
#include <cfloat>
#include <iomanip>

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class Clone_gof { 
private:
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	float_t number_;
	std::string name_;

	Clone_gof();	
	Clone_gof(std::vector <std::string>) : Clone_gof(){};	
	Clone_gof(const std::string &, const float_t &);	

	std::string header(void) const;
	size_t size(void) const;

	friend std::ostream& operator<< (std::ostream&, const Clone_gof&);	//!< use the << operator to write allele_stat.
	friend std::istream& operator>> (std::istream&, Clone_gof&);		//!< use the >> operator to read allele_stat.
	static const std::string file_name;					//!< The dafualt extention for files.
	static const std::string table_name;					//!< Destination table in Db.
};

#endif
