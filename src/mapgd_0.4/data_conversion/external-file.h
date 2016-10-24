#ifndef _EXTERNAL_FILE_H_
#define _EXTERNAL_FILE_H_	

#include <map>
#include <string>
#include <iostream>
#include <typeinfo>

#include "stream-tools.h"
#include "typedef.h"
#include "stream-tools.h"
#include "map-file.h"

/* These should all be included from the objects using this code, but I'm not 
 * linking things correctly right now.
 */
#include "vcf-file.h"

// PLEASE LIMIT LINE LENGTH TO 79 CHARACTERS----------------------------------/

//! An interface for reading and writing data specified outside of mapgd.
/*! Because we do not know the interanl structure of the external data, we 
 * ...
 * 
 * in the table, which must store a pair of IDs retrived by the get_id0() and 
 * get_id1() member functions. Additionally, data types must return a name for 
 * the table in the database where the data may be stored with member function 
 * table_name(), and must name the columns of the table with cannonical data 
 * types listed in the file ?. This allows data to easily be transfered into 
 * and out of a database. Finally...
 *
 */

template <class External_data>
class External_file : public Base_file {
private :
protected :
	void read_binary( External_data &);		//!< Read file in binary mode.		
	void write_binary(const External_data &);	//!< Write in binary mode.

	virtual void read_text( External_data &){};	//!< Read file in text mode.
	virtual void write_text(const External_data&){};//!< Write file in text mode.

	using Base_file::out_;
	using Base_file::in_;

	using Base_file::open_;
	using Base_file::read_;
	using Base_file::write_;
	using Base_file::binary_;
	using Base_file::filename_;
public:
	using Base_file::open;

	///! The function that opens a External_file (if file).
	void open(const char *, const std::ios_base::openmode &);

	//! Writes a row to the file and advances one row. 
	/*! Returns the ostream.
	 */
	External_file& write(const External_data &);

	//! Reads a row from the file and advances one row. 
	/*! Returns the class.
	 */
	External_file& read(External_data &);				

//	void close(void);
};

#endif
