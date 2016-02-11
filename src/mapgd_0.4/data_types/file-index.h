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
#include "data.h"

/// An interface that transforms pairs of name and position keys into record numbers.
/* File_index should always be used to re designate pairs of keys into a single record
 * number. 
 */
class File_index : public virtual Data {

private:
	//! a hash table to store pairs of strings and the numeral associated with them.
	std::map <std::string, id0_t> id0_str_;			

	//! the keys to the hash table id0_str_.
	std::vector <std::string> id0_;	
			
	//! the number of bp per scaffold.
	std::vector <id1_t> size_;				

	//! the cumulative number of bp per scaffold.
	std::vector <id1_t> cumulative_size_;	

	//! the last numberical id returned by an decodeid0 querry. Initilized to 0-1.
	id0_t last_id0_;						

	//! the last string id returned by an decodeid0 querry. Initilized to "".
	std::string last_id0_str_;
	
	//! TODO: Is this ever used?
	bool open_;

	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new File_index(Columns);
	}

public:
	File_index();						
	File_index(std::vector<std::string>) : File_index(){};						

	//! returns the number of rows until a row with id0, id1.
	id1_off_t get_offset (const std::string &, const id1_t &) const;	

	//! returns the record number corresponding to id0, id1.
	id1_t get_abs_pos (const std::string &, const id1_t &) const;	

	//! returns the record number corresponding to id0, id1.
	id1_t get_abs_pos (const id0_t &, const id1_t &) const;	

	//! returns the id0 corresponding to string.
	id0_t get_id0 (const std::string &) const;	

	//! returns the id0 corresponding to rowid.
	id0_t get_id0 (const id1_t &) const;	

	//! returns the id1 corresponding rowid.
	id1_t get_id1 (const id1_t &) const;	

	//! returns the string corresponding to id0.
	std::string get_string (const id0_t &) const;	

	//! sets the string corresponding to id0.
	void set_string (const id0_t &, const std::string &);	

	//! set the number of rows (bp) in the scaffold with string id id0.
	void set_size (const std::string &, const id1_t &);	

	//! set the number of rows (bp) in the scaffold with numeric id id0.
	void set_size (const id0_t &, const id1_t &);	

	//! set the number of rows in the most resently read scaffold.
	void set_last_size (const id1_t &);			

	//! set the number of rows in the next scaffold to be read.
	void set_next_size (const id1_t &);			

	//! get a vector containing the sizes of each scaffold.
	std::vector <id1_t> get_sizes (void);		

	//! get the size of scaffold blarg.
	id1_t get_size (const id0_t &) const;

	/// get the size of scaffold blarg.
	id1_t get_size (const std::string &) const;

	//! get the sum of all size less than or equal to id0.
	id1_t get_cumulative_size (const id0_t &) const;

	//! get the sum of all size less than or equal to id0.
	id1_t get_cumulative_size (const std::string &) const;

	//! get the sum of all size.
	id1_t get_reference_size (void) const;

	//! get the size of the last scaffold which has been read from. returns map_file::noid on error.
	id1_t get_last_size (void);				

	//! get the size of the next scaffold to be read. returns map_file::noid on error.
	id1_t get_next_size (void);				

	//! read in a File_index (strictly text mode).
	int read_index(std::istream &);				

	//! read in a file index from a sam header.
	std::istream & from_sam_header(std::istream &);

	//! write out a File_index
	int write_index(std::ostream &);			

	//! not implemented. Do not use.
	File_index & operator=(const File_index&);

	bool is_open(void) const;

	std::string header(void) const;
	static const std::string file_name;
	static const std::string table_name;
	size_t size(void) const;

	//! use the << operator to write File_index.
	void write (std::ostream&) const;	
	//! use the >> operator to read File_index.
	void read (std::istream&);	
};

#endif 
