#ifndef FILE_INDEX_H_
#define FILE_INDEX_H_

#include "stream_tools.h"
#include "typedef.h"
#include <iostream>
#include <map>
#include <vector>
#include <stdio.h>
#include <cstring>

/// An interface that transforms pairs of name and position keys into record numbers.
/* file_index should always be used to re designate pairs of keys into a single record
 * number. 
 */
class file_index{

private:
/*! \breif A hash table to store pairs of strings and the numeral associated with them.
 */
	std::map <std::string, id0_t> id0_str_;			
/*! \breif The keys to the hash table id0_str_.
 */
	std::vector <std::string> id0_;	
			
/*! \breif The sizes of the rows in bytes.
 */
	std::vector <id1_t> size_;				

/*! \breif The last numberical id returned by an decodeid0 querry. Initilized to 0-1.
 */
	id0_t last_id0_;						

/*! \breif The last string id returned by an decodeid0 querry. Initilized to "".
 */
	std::string last_id0_str_;
	char delim_;				
	char key_;			

public:
	file_index();						

	/// returns the number of rows until a row with id0, id1.
	id1_t get_offset (const std::string &, const id1_t &);	

	/// rerutnrs the record number corresponding to id0, id1.
	id1_t get_row (const std::string &, const id1_t &);	

	/// set the number of rows (bp) in the scaffold with string id id0.
	void set_size (const std::string &, const id1_t &);	

	/// set the number of rows (bp) in the scaffold with numeric id id0.
	void set_size (const id0_t &, const id1_t &);	

/*! \breif Set the number of rows in the most resently read scaffold.
 */
	void set_last_size (const id1_t &);			

/*! \breif Set the number of rows in the next scaffold to be read.
 */
	void set_next_size (const id1_t &);			

/*! \breif Get a vector containing the sizes of each scaffold.
 */
	std::vector <id1_t> get_all_sizes (void);		

/*! \breif Get the size of the last scaffold read.
 */
	id1_t get_last_size (void);				

/*! \breif Get the size of the next scaffold to be read.
 */
	id1_t get_next_size (void);				

/*! \breif Read in a file_index (strictly text mode).
 */
	int read_index(std::istream &);				

/*! \breif Read in a file_index from a sam header.
 */
	int from_sam_header(std::istream &);

/*! \breif Write out a file_index (strictely text mode).
 */
	int write_index(std::ostream &);			

/*! \breif Not implemented. Do not use.
 */
	file_index & operator=(const file_index&);
};

#endif 
