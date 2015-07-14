//TODO, these classes need to be split apart...
#ifndef GCF_FILE_H_
#define GCF_FILE_H_	
#include <map>
#include <string>

#include "typedef.h"
#include "stream-tools.h"
#include "genotype.h"
#include "gcf-file-header.h"

/* \breif used to index and deindex files*/

	
class  gcf_file {
//private variables should be initialized by reading the header...
private:
	bool open_;					// indicates whether the profile opened succesfully

	/*these should all be controlled through the header_*/

	char delim_column;				// the delimiter which seperates columns
	char delim_genotype;				// the delimiter which seperates genotypes

	unsigned int samples_;				// the number of samples (i.e. different individuals or populations) in the profile.
	count_t size_;					// the number of lines in the sample. 0 if unkown.

	bool read_;					// file is open for reading.
	bool write_;					// file is open for writing.
	bool binary_;					// binary mode flag. Incompatable with mpileup and noheader flags.

	/*done*/

	static const std::string names_;		// ACGTN

	int readt(int);					//read file in text mode.
	int readb(int);					//read file in binary mode.

	int writet();					//write stat information in memory to file in text mode.
	int writeb();					//write stat information in memory to file in binary mode.

	file_index *index_;

	std::istream *in;				// all data is read from in.
	std::ostream *out;				// all data is writen is writen to out.
	std::fstream inFile;				// the file to read data from (if not stdin).
	std::ofstream outFile;				// the file to write data to (if not stdout).
public:
	gcf_file();					//default constructor

	gcf_file* open(const char *, const char *);	//The function that opens a profile (if file).
	gcf_file* open(const char *);			//The function that opens a profile (if stdin).
	bool is_open(void) const;			//returns true if profile is open, false otherwise.

	/*basic io operation*/
	int copy(const gcf_file&);			//copys a line from profile
	int read();					//reads a line from the instream. Returns 0 on success, EOF on EOF.
	int write();					//writes a line to the outstream. Returns 0 on success, EOF on EOF.

	void close(void);				//close iostreams, writes tail, etc.

	/*functions dealing with the header*/

	void set_delim_column(const char&);		// the delimiter which seperates columns
	void set_delim_genotype(const char&);		// the delimiter which seperates columns

	int set_samples(count_t);			//set the number of samples in the profile (only called in write mode).
	int set_columns(count_t);			//set the number of columns for reading and writing.

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
	count_t size(void) const;			//number of populations/individuals

	const count_t get_index(count_t) const;		//returns the index of the alleles in order a sorted order

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
