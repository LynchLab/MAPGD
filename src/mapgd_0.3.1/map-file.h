/* TODO: Implement.
 * 
 */
#ifndef MAP_FILE_H_
#define MAP_FILE_H_	

#include <map>
#include <string>

#include "typedef.h"
#include "stream-tools.h"
#include "genotype.h"
#include "file-index.h"
#include "map-file-header.h"

/// The file used to pass data between mapgd commands. 
/* A map file is a table with columns specified by the data types stored in the table, and row specified by a row number. 
 * It has only read/write possition, which is the row number of the row being manipulated. This may be replaced by some form of
 * MySQL request/strucutre/something, at which point this will just be a wrapper.
 *
 */
class  map_file {
//private variables should be initialized by reading the header...
private:
	bool open_;					//!< indicates whether the profile opened succesfully

	/*these should all be controlled through the header_*/

	char delim_column_;				//!< The delimiter which seperates columns
	bool read_;					//!< File is open for reading.
	bool write_;					//!< File is open for writing.
	bool binary_;					//!< Binary mode flag. Incompatable with mpileup and noheader flags.

	/*done*/

	std::istream * pin_;				//!< All data is read from in.
	std::ostream * pout_;				//!< All data is writen is writen to out.

	std::fstream file_;				//!< The file for reading/writing data (if not stdin).

	int readt();					//!< Read file in text mode.
	int readb();					//!< Read file in binary mode.

	int writet();					//!< Write stat information in memory to file in text mode.
	int writeb();					//!< Write stat information in memory to file in binary mode.
public:
	map_file();					//!< default constructor

	bool is_open(void) const;			//!< Returns true if profile is open, false otherwise.

	/** @defgroup BasicIO Basic IO Operations
	 * @{
  	 */
	void map_seekg(id1_t pos);				//!< Seeks to row number for gets. 
	void map_seekp(id1_t pos);				//!< Seeks to row number for puts.

	void map_seekg(id1_t off, std::ios_base::seekdir way);	//!< Seeks to row number for puts.
	void map_seekp(id1_t off, std::ios_base::seekdir way);	//!< Seeks to row number for gets. 

	id1_t tellp(void);					//!< Tells row number of puts.
	id1_t tellg(void);					//!< Tells row number of gets.

	map_file* open(const char *, const char *);		//!< The function that opens a mapfile (if file). These should return void, but I just can't help myself.
	map_file* open(const char *);				//!< The function that opens a mapfile (if stdin). These should return void, but I just can't help myself.
	void close(void);					//!< Close iostreams, writes tail, etc.

	/** @} */

	/*functions dealing with the header*/

	/*functions dealing with ?*/
	size_t size(void) const;				//!< Retuern the number of rows in the file. Returns 0 if unknown.
								// If pro file is opened from a std::in, then the number of rows is simply 
								// the number of rows that have been read in until d. However, if 

	size_t column_size(void) const;				//!< Returns the number of columns in the file. Returns 0 if unknown.

	file_index get_index(void) const;		//!< Returns the file_index.
	const id1_t get_line_number(void) const;	//!< 

};
	
map_file & read_row(id0_t, map_file, );			//!< Reads a row from the file and advances the read one row. Returns 0 on success, EOF on EOF.

map_file & write_row(,);				//!< Writes a row to the file and advances one row. Returns 0 on success, EOF on EOF.
							//
	
#endif
