#ifndef FILE_INDEX_H_
#define FILE_INDEX_H_

#include "file-index.h"

/*! \breif and index file designed to help with serializiation and 
 */

class file_index{

private:
	std::map <std::string, count_t> id0_str_;		//
	std::vector <std::string> id0_;				//
	std::vector <id1_t> sizes_;				//

	id0_t lastid0_;						//initilize to 0-1;
	std::string lastid0_str_;				//initilize to "";

public:
	file_index();						// does not initilize . . .

	const id0_t & encode_id0(const std::string &);		//
	const std::string & decode_id0(const id0_t &);		//

	id1_t get_offset (const id0_t &, const id1_t &);	//

	void set_size (const std::string &, const id1_t &);	//
	void set_last_size (const id1_t &);			//
	void set_next_size (const id1_t &);			//

	std::vector <id1_t> get_all_sizes (void);		//

	id1_t get_last_size (void);				//
	id1_t get_next_size (void);				//

	int read_index(std::istream *);				// . . .
	int write_index(std::ostream *);			// . . .

	file_index & operator=(const file_index&); 		//I don't think the copy makes sense. . . 
};

#endif 
