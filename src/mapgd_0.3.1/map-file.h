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

//ios::beg	offset counted from the beginning of the stream
//ios::cur	offset counted from the current position
//ios::end	offset counted from the end of the stream


/// The file used to pass data between mapgd commands. 
/* A map file is a table with colulmns specified by the data types stored in the table, and row specified by a row number. 
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

	row_index *row_index_;				//!< An row_index which transforms to locate rows in the file.

	std::istream *in_;				//!< All data is read from in.
	std::ostream *out_;				//!< All data is writen is writen to out.

	std::fstream in_file_;				//!< The file to read data from (if not stdin).
	std::ofstream out_file_;				//!< The file to write data to (if not stdout).

	int readt();					//!< Read file in text mode.
	int readb();					//!< Read file in binary mode.

	int writet();					//!< Write stat information in memory to file in text mode.
	int writeb();					//!< Write stat information in memory to file in binary mode.
public:
	map_file();					//!< default constructor

	map_file* open(const char *, const char *);	//!< The function that opens a mapfile (if file).
	map_file* open(const char *);			//!< The function that opens a mapfile (if stdin).
	void close(void);				//!< Close iostreams, writes tail, etc.

	bool is_open(void) const;			//!< Returns true if profile is open, false otherwise.

	/** @defgroup BasicIO Basic IO Operations
	 * @{
  	 */

	seekp(id1_t);					//!< Seeks to row number 
	seekg(id1_t);					//!< Seeks to row number 
	tellg();					//!< 
	tellp();					//!< 

	/** @} */

	/*functions dealing with the header*/

	void set_delim_column(const char&);			//!< Sets the delimiter that seperates columns. Only used in text mode.
	const char & get_delim_column(const char&) const;	//!< Gest the delimiter that seperates columns. Only used in text mode.

	/*functions dealing with ?*/
	size_t size(void) const;			//!< Retuern the number of rows in the file. Returns 0 if unknown.
	size_t column_size(void) const;			//!< Returns the number of columns in the file. Returns 0 if unknown.

	file_index get_index(void) const;		//!< Returns the file_index.
	const count_t get_line_number(void) const;

};
	
map_file & read(id0_t, map_file, );			//!< Reads a row from the file and advances the read one row. Returns 0 on success, EOF on EOF.

map_file & write(,);					//!< Writes a row to the file and advances one row. Returns 0 on success, EOF on EOF.
							//
	
#endif
