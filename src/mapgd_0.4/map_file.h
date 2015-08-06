/* TODO: Implement.
 * 
 */
#ifndef MAP_FILE_H_
#define MAP_FILE_H_	

#include <map>
#include <string>

#include "datatypes.h"
#include "typedef.h"
#include "stream_tools.h"
#include "file_index.h"

/// The pipe used to pass data between mapgd commands. 
/* It is a purely serial interface (i.e. you cannot use seekg or seekp). A tiny bit of information can
 * be passed through the keys, but keys can only be . . . 
 */
class  map_pipe {
//private variables should be initialized by reading the header...
private:
	bool open_;					//!< indicates whether the profile opened succesfully

	/*these should all be controlled through the header_*/

	char delim_column_;				//!< the delimiter which seperates columns
	bool read_;					//!< pipe is open for reading.
	bool write_;					//!< pipe is open for writing.
	bool binary_;					//!< binary mode flag. Incompatable with mpileup and noheader flags.

	/*done*/

	std::istream * in_;				//!< all data is read from in.
	std::ostream * out_;				//!< all data is writen is writen to out.

	std::fstream file_;				//!< the file for reading/writing data (if not stdin).

public:
	map_pipe();					//!< default constructor

	bool is_open(void) const;			//!< returns true if profile is open, false otherwise.

	/** @defgroup BasicIO Basic IO Operations
	 * @{
  	 */

	map_pipe* open(const char *, std::ios_base::openmode);	//!< the function that opens a mapfile (if file).
	map_pipe* open(std::ios_base::openmode);		//!< the function that opens a mapfile (if stdin).
	void close(void);					//!< close iostreams.

	/** @} */

	/** @defgroup BasicKey key
	 * @{
  	 */

	std::list <key> get_keys(void);			//!< returns all the keys in the mapfile.

	void format(std::list <key>);			//!< formats the file. 
							//This can only be done if row is not initalized.

	/** @} */
	size_t width(void) const;			//!< Returns the number of bytes in a row.

	file_index get_index(void) const;		//!< Returns the file_index.
};
	
int read_row(map_pipe &, row &);		//!< Reads a row from the file and advances the read one row. 
						// Returns 0 on success, EOF on EOF.

int write_row(map_pipe &, const row &);		//!< Writes a row to the file and advances one row. 
						//Returns 0 on success, EOF on EOF.
	
#endif
