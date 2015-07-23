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

/*! \breif used for passing information around within the program.
 *
 */
class  map_file {
//private variables should be initialized by reading the header...
private:
	bool open_;					//!< indicates whether the profile opened succesfully

	/*these should all be controlled through the header_*/

	char delim_column;				//!< The delimiter which seperates columns
	char delim_genotype;				//!< The delimiter which seperates genotypes

	unsigned int samples_;				//!< The number of samples (i.e. different individuals or populations) in the profile.
	count_t size_;					//!< The number of lines in the sample. 0 if unkown.

	bool read_;					//!< File is open for reading.
	bool write_;					//!< File is open for writing.
	bool binary_;					//!< Binary mode flag. Incompatable with mpileup and noheader flags.

	/*done*/

	static const std::string names_;		//!< ACGTN

	int readt(int);					//!< Read file in text mode.
	int readb(int);					//!< Read file in binary mode.

	int writet();					//!< Write stat information in memory to file in text mode.
	int writeb();					//!< Write stat information in memory to file in binary mode.

	file_index *index_;

	std::istream *in;				//!< All data is read from in.
	std::ostream *out;				//!< All data is writen is writen to out.
	std::fstream inFile;				//!< The file to read data from (if not stdin).
	std::ofstream outFile;				//!< The file to write data to (if not stdout).
public:
	gcf_file();					//!< default constructor

	gcf_file* open(const char *, const char *);	//!< The function that opens a mapfile (if file).
	gcf_file* open(const char *);			//!< The function that opens a mapfile (if stdin).
	bool is_open(void) const;			//!< Returns true if profile is open, false otherwise.

	/*basic io operation*/
	int copy(const gcf_file&);			//!< Copys a line from profile
	int read();					//!< Reads a line from the instream. Returns 0 on success, EOF on EOF.
	int write();					//!< Writes a line to the outstream. Returns 0 on success, EOF on EOF.

	void close(void);				//!< Close iostreams, writes tail, etc.

	/*functions dealing with the header*/

	void set_delim_column(const char&);		//!< The delimiter which seperates columns
	void set_delim_genotype(const char&);		//!< The delimiter which seperates columns

	int set_samples(count_t);			//!< Set the number of samples in the profile (only called in write mode).
	int set_columns(count_t);			//!< Set the number of columns for reading and writing.

	int set_column_name(const count_t&, const std::string &);
	int set_sample_name(const count_t&, const std::string &);

	const std::string get_sample_name(const count_t &) const;

	/*a set of functions for converting the string information read in text mode into count_t and vice versa.*/

	const id0_t encode_id0(const std::string &);
	const id1_t encode_id1(const std::string &);

	const std::string decode_id0(const count_t &);
	const std::string decode_id1(const uint64_t &);
	const std::string decode_extraid(const count_t &, const count_t &);

	/*functions dealing with ?*/
	count_t size(void) const;			//!< Number of populations/individuals

	const count_t get_index(count_t) const;		//!< Returns the index of the alleles in order a sorted order

	std::string get_ids(void); 			//

	const count_t get_line_number(void) const;

	const count_t get_id0(void) const;
	const uint64_t get_id1(void) const;
	const count_t get_extraid(const count_t &) const;

	void set_id0(const count_t &);
	void set_id1(const uint64_t &);
	void set_extraid(const count_t &, const count_t &);
};
	
#endif
