/* synonym for population? */

#ifndef SAMPLE_GOF_H
#define SAMPLE_GOF_H

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "typedef.h"
#include "data.h"
#include "stream_tools.h"

///	The sample specific goodness of fit.
class Sample_gof : public Data { 
private:
	void write (std::ostream&) const;	//!< use the << operator to write Allele.
	void read (std::istream&);		//!< use the >> operator to read Allele.
	void sql_read (std::istream&);		//!< use the >> operator to read Allele.
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Sample_gof(Columns);
	}
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	float_t number_;
	count_t smp_num_;
	std::string name_;
	

	Sample_gof();	
	Sample_gof(std::vector <std::string>) : Sample_gof(){};	
	Sample_gof(const std::string &, const float_t &);	
	Sample_gof(const count_t &, const std::string &, const float_t &);	

	std::string header(void) const;
	size_t size(void) const;

	static const std::string file_name;				//!< The dafualt extention for files.
	static const std::string table_name;				//!< Destination table in Db.

	const std::string get_file_name(void) const;			//!< The dafualt extention for files.
	const std::string get_table_name(void) const;			//!< Destination table in Db.

	static const bool binary;					//!< Destination table in Db.

	const bool get_binary () const;					//!< Destination table in Db.

	const std::string sql_header(void) const;	//!< Returns a string to create an SQL table.
	const std::string sql_column_names(void) const; //!< Returns the column names in the SQL table.
	const std::string sql_values(void) const;	//!< Returns a string to insert values in an SQL table.
};

#endif
