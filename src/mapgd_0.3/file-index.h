#ifndef _FILE_INDEX_H_
#define _FILE_INDEX_H_


//these need to be moved ...

#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>	
#include <cstring>	//memcpy
#include <string>

#include "stream-tools.h"
#include "typedef.h"


//#include "datatypes/row.h"

/// An interface that transforms pairs of name and position keys into record numbers.
/* file_index should always be used to re designate pairs of keys into a single record
 * number. 
 */
class file_index{

private:
	/// a hash table to store pairs of strings and the numeral associated with them.
	std::map <std::string, id0_t> id0_str_;			

	/// the keys to the hash table id0_str_.
	std::vector <std::string> id0_;	
			
	/// the sizes of the rows in bytes.
	std::vector <id1_t> size_;				

	/// the last numberical id returned by an decodeid0 querry. Initilized to 0-1.
	id0_t last_id0_;						

	/// the last string id returned by an decodeid0 querry. Initilized to "".
	std::string last_id0_str_;
	
	bool open_;

	/// this is the stupid implementation for now, make it better later.
	//gap gaps_; //You know what? F$%# it. It will make our files 20% larger, do I really care that much?

public:
	file_index();						

	/// returns the number of rows until a row with id0, id1.
	id1_off_t get_offset (const std::string &, const id1_t &) const;	

	/// returns the record number corresponding to id0, id1.
	id1_off_t get_rowid (const std::string &, const id1_t &) const;	

	/// returns the id0 corresponding to string.
	id0_t get_id0 (const std::string &) const;	

	/// returns the id0 corresponding to rowid.
	id0_t get_id0 (const id1_t &) const;	

	/// returns the id1 corresponding rowid.
	id1_t get_id1 (const id1_t &) const;	

	/// returns the string corresponding to id0.
	std::string get_string (const id0_t &) const;	

	/// sets the string corresponding to id0.
	void set_string (const id0_t &, const std::string &);	

	/// set the number of rows (bp) in the scaffold with string id id0.
	void set_size (const std::string &, const id1_t &);	

	/// set the number of rows (bp) in the scaffold with numeric id id0.
	void set_size (const id0_t &, const id1_t &);	

	/// set the number of rows in the most resently read scaffold.
	void set_last_size (const id1_t &);			

	/// set the number of rows in the next scaffold to be read.
	void set_next_size (const id1_t &);			

	/// get a vector containing the sizes of each scaffold.
	std::vector <id1_t> get_sizes (void);		

	/// get the size of scaffold blarg.
	id1_t get_size (const id0_t &) const;

	/// get the size of the last scaffold which has been read from. returns map_file::noid on error.
	id1_t get_last_size (void);				

	/// get the size of the next scaffold to be read. returns map_file::noid on error.
	id1_t get_next_size (void);				

	/// read in a file_index (strictly text mode).
	int read_index(std::istream &);				

	/// read in a file index from a sam header.
	int from_sam_header(std::istream &);

	/// write out a file_index
	int write_index(std::ostream &);			

	/// not implemented. Do not use.
	file_index & operator=(const file_index&);

	bool is_open(void) const;
};

#endif 
