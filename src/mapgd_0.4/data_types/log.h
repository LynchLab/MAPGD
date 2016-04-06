/* synonym for population? */

#ifndef _LOG_H_
#define _LOG_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>

#include "typedef.h"
#include "data.h"
#include "stream-tools.h"

///	A class to store a log of commands that have been run.
/** Commands should set logs whenever they are run.
 */
class Log : public virtual Data { 
private:
	void read(std::istream& str);
	///! the write function must be ? by the child class.
	void write(std::ostream& str) const;
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Log(Columns);
	};
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	std::string command;
	std::string message;
	time_t time;

	Log();	
	Log(const std::vector <std::string> &) : Log(){};	
	Log(const std::string &, const float_t &);	

	static const std::string file_name;				//!< The dafualt extention for files.
	static const std::string table_name;				//!< Destination table in Db.

	const std::string header(void) const;
	size_t size(void) const;

	const std::string get_file_name(void) const;				//!< The dafualt extention for files.
	const std::string get_table_name(void) const;				//!< The dafualt extention for files.
	const std::string sql_header(void) const;				//!< Destination table in Db.
	const std::string sql_column_names(void) const;				//!< Destination table in Db.
	const std::string sql_values(void) const;				//!< Destination table in Db.
};

//! Get current date/time, format is YYYY-MM-DD.HH:mm:ss, by TrungTN
const std::string format_time(cosnt time_t &);
#endif
