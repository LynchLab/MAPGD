/* synonym for population? */

#ifndef _SAMPLE_NAME_H_
#define _SAMPLE_NAME_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "data.h"
#include "stream_tools.h"

///	The name of a sample (i.e. what the user called it).
class Sample_name : public virtual Data { 
private:
	void read(std::istream& str);
	///! the write function must be ? by the child class.
	void write(std::ostream& str) const;
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Sample_name(Columns);
	};
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	std::string mpileup_name;
	std::vector <std::string> sample_names;

	Sample_name();	
	Sample_name(const std::vector <std::string> &) : Sample_name(){};	
	Sample_name(const std::string &, const float_t &);	

	static const std::string file_name;				//!< The dafualt extention for files.
	static const std::string table_name;				//!< Destination table in Db.
	static const bool binary;				//!< Destination table in Db.

	std::string header(void) const;
	size_t size(void) const;

	const std::string get_file_name(void) const;				//!< The dafualt extention for files.
	const std::string get_table_name(void) const;				//!< The dafualt extention for files.
	const bool get_binary(void) const;					//!< The dafualt extention for files.

	const std::string sql_header(void) const;				//!< Destination table in Db.
	const std::string sql_column_names(void) const;				//!< Destination table in Db.
	const std::string sql_values(void) const;				//!< Destination table in Db.
};

#endif
