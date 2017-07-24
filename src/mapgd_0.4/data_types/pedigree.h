#ifndef _PEDIGREE_DATA_H_
#define _PEDIGREE_DATA_H_

#include <cstring>

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "typedef.h"
#include "data.h"

#define E_LIM 25

/// Pedigree data.
class Pedigree : public Data{ 
private:
	void write (std::ostream&) const;	//!< use to write Allele. Inherits <<
	void read (std::istream&);		//!< use to read Allele. Inherits >>
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Pedigree(Columns);
	}
public:
	Pedigree_record *parents;
	Pedigree_record *roots;

	char delim;	//!< the delimiter used when reading/writing the class in text mode.	

	Pedigree();	
	//! Delegating a neccisary constructor.	
	Pedigree(const std::vector <std::string> &) : Pedigree(){}; 
	//! Construct with names. 
	Pedigree(const std::string &, const std::string &);		  

	//! The header line of plain text files. 
	std::string header(void) const;
	//! Size in bytes for binary read/write. 
	size_t size(void) const;

	static const std::string file_name;	//!< The dafualt extention for files.
	static const std::string table_name;	//!< Destination table in Db.
	static const bool binary;

	const bool get_binary() const;

	Pedigree& operator=(const Pedigree &rhs);
};

#endif
