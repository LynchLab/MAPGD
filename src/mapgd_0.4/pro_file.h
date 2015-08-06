#ifndef PROFILE_H_
#define PROFILE_H_	

#include <cstdio>
#include <cstring>
#include <iostream>
#include <climits>
#include <stdlib.h> 	//Needed for one call to atoi

#include <list>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <ios>

#include "typedef.h"
#include "stream_tools.h"
#include "file_index.h"
#include "datatypes.h"

#define NONE		0
#define BADHEADER	1
#define UNEXPECTED	2
#define VER		2.0

#define SKIP		1

#define NEWID0		1
#define EOBIN		2

///breif A table to store quartet information. 
/* 
 *	Eventually this will be a specialization to add conversion functions from mpileup, and the various .pro 
 *	format to a more map format that will be used to read/write information within mapgd commands. 
 *	See the design document for the map format in the docs directory. 
 */
class  profile {
private:
	/** @defgroup pro_format .pro format options.
	 * @{
	 */
	char delim_column;			//!< the delimiter which seperates columns
	char delim_quartet;			//!< the delimiter that seperates counts in a quartet
	size_t columns_;			//!< 5|6|7|more?
	
	/** @} */

	size_t samples_;			//!< the number of samples (i.e. different individuals or populations) in the profile.
	id1_t size_;				//!< the number of lines in the sample. 0 if unkown.

	bool open_;				//!< indicates whether the profile opened succesfully
	bool sig_;				//!< were alleles thrown out if the allele only occurred in reads from one direction?
	bool read_;				//!< file is open for reading.
	bool write_;				//!< file is open for writing.
	bool binary_;				//!< binary mode flag. Incompatable with mpileup and noheader flags.
	bool mpileup_;				//!< file is an mpileup. Setting this should set the noheader flag.
	bool noheader_;				//!< file has no header

	locus *locus_;				//!< a structure that stores quartet information.
	gt_t  *ref_;				//!<
	id1_t *rowid_;				//!<

	id0_t id0_, last_id0_;			//!< something is needed to support Berhard's format.
	id1_t id1_;

	int readm(row &);			//!< read file in mpileup mode.
	int readt(row &);			//!< read file in text mode.

	int writet();				//!< write the quartet information in memory to file in text mode.
	int writet(row &);		//!< write the quartet information passed to file in text mode.

	file_index index_;			//!< this tells us where we are.

	std::istream *in;			//!< all data is read from in.
	std::ostream *out;			//!< all data is writen is writen to out.
	std::fstream file;			//!< the file to read/write data to/from (if not stdin/stdout).

	int readheader(void);			//!< attempts to format the profile correctly. 
	void inline scan(const std::string &, quartet_t &); 
public:
	profile();					//!< default constructor
	static constexpr char *names_="ACGTN";		//!< ACGTN
	void open(const char *, std::ios_base::openmode);	//!< The function that opens a profile (if file).
	void open(std::ios_base::openmode);			//!< The function that opens a profile (if stdin).

	bool is_open(void) const;			//!< Returns true if profile is open, false otherwise.

	/*basic io operation*/
	int read(row&);					//!< reads a row (i.e. line) from the istream associated w/ profile. 
							//* Returns 0 on successs, EOF on EOF. Exist if no stream is open. */

	int write(row&);				//!< writes a row (i.e. line) to the ostream associated w/ profile.
							//* Returns 0 on successs, EOF on EOF. Exist if no stream is open. */

	void close(void);				//!< close iostreams, writes tail, etc.

	void static merge(std::list <profile *>);	//!< take multiple profile and combine them into a single profile.

	/*functions dealing with the header*/
	/** @defgroup pro_format .pro format options.
          */
	void set_delim_column(const char&);		//!< The delimiter which seperates columns
	void set_delim_quartet(const char&);		//!< The delimiter which seperates columns

	void set_samples(const size_t &);			//!< Set the number of samples in the profile (only called in write mode).
	void set_columns(const size_t &);			//!< Set the number of columns for reading and writing.
	size_t get_columns(void) const {return columns_;};	//!< Set the number of columns for reading and writing.

	int set_column_name(const size_t&, const std::string &);
	int set_sample_name(const size_t&, const std::string &);
	int set_sample_name(const std::string &);

	const std::string get_sample_name(const size_t &) const;

	/*functions dealing with ?*/
	size_t size(void) const;			//!< Returns the number of samples/individuals in a profile.
	size_t length(void) const {return size_;};	//!< Returns the number of rows in a profile.

	/*functions dealing with the quartets*/

	void set_base_count(const size_t &x, const gt_t &y, const count_t &z);	//!< Sets sample x, base y to z.

	/*playing with the mask*/

	void mpileup(void) {mpileup_=true; noheader_=true;}	//!< File is an mpileup. Setting this should set the noheader flag.
	void noheader(void) {noheader_=true;}			//!< File has no header

	size_t masked_count(void);			//!< Returns a count of the number of quartets masked at the current locus. 
	void maskall(void);				//!< Mask all lines at the current locus. 


	void header(void) {noheader_=false;};		//!< If header is true pro-file will expect a file to have a header...

	//locus

	locus get_locus(void) {return *locus_;};		//!< Returns the current locus.
	void set_locus(const locus &locus) {*locus_=locus;};	//!< Set the current locus to locus.
	void sort(void);				//!< Left sort the sum of unmaksed quartets (i.e. record 0 is the larges).

	//quartet

	quartet_t get_quartet(const size_t &t) {return locus_->get_quartet(t);};	//!< Returns the quartet of sample t at the current locus.
	std::vector <quartet_t>::iterator begin(void) {return locus_->begin();};	//!< Returns an iterator to the quartet at the begining of the current locus.
	std::vector <quartet_t>::iterator end(void) {return locus_->end();};	//!< Returns an iterator to the quartet just past the end of the current locus. */

};
	
	
#endif
