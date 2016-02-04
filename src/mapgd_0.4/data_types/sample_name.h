/* synonym for population? */

#ifndef _SAMPLE_NAME_H_
#define _SAMPLE_NAME_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "../typedef.h"
#include "data.h"
#include "../stream-tools.h"

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class Sample_name : public virtual Data { 
private:
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	std::string mpileup_name;
	std::vector <std::string> sample_names;

	Sample_name();	
	Sample_name(std::vector <std::string>) : Sample_name(){};	
	Sample_name(const std::string &, const float_t &);	

	static const std::string file_name;				//!< The dafualt extention for files.
	static const std::string table_name;				//!< Destination table in Db.

	std::string header(void) const;
	size_t size(void) const;

	friend std::ostream& operator<< (std::ostream&, const Sample_name&);	//!< use the << operator to write allele_stat.
	friend std::istream& operator>> (std::istream&, Sample_name&);		//!< use the >> operator to read allele_stat.

	const std::string get_file_name(void) const;				//!< The dafualt extention for files.
	const std::string get_table_name(void) const;				//!< The dafualt extention for files.
	const std::string sql_header(void) const;				//!< Destination table in Db.
	const std::string sql_column_names(void) const;				//!< Destination table in Db.
	const std::string sql_values(void) const;				//!< Destination table in Db.
};

#endif
