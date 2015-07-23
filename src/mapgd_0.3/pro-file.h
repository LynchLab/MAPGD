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

#include "typedef.h"
#include "stream-tools.h"
#include "locus.h"
#include "quartet.h"

#define NONE		0
#define BADHEADER	1
#define UNEXPECTED	2
#define VER		2.0

#define SKIP		1

#define NEWID0		1
#define EOBIN		2

class profile;

struct profile_flags{
};


/*! \breif A class to control the varrious atributes of .pro files.  A class to control the various attributes of .pro files. This class was a terrible design decision that I am now paying for. Some form of this class will likely persist, but in the future I think it will only set internal flags. 
 */
class  profile_header{
private:
//private variables should be initialized by reading the header...
	std::map <std::string, count_t> id0_str;
	std::vector <std::string> id0;
	std::vector <std::string> extraids;		//!< extra ids associated with the quartet. (ref base identiy?).	

	std::vector <std::string> column_names;		//!< the names of all the columns in the profile.
	
	count_t lastid0;				//!< initilize to 0-1;
	std::string lastid0_str;			//!< initilize to "";

	count_t encodechar[256];
	char decodechar[256];

	bool *sig_;					// were alleles thrown out if the allele only occurred in reads from one direction?
	bool *read_;					// ?
	bool *write_;					// ?
	bool *binary_;					// ?
	bool *mpileup_;					// ?
	bool *noheader_;				// ?

	id1_t *size_;					// the number of lines in the sample. 0 if unkown.

	char *delim_column;				// the delimiter that seperates columns.
	char *delim_quartet;				// the delimiter that seperates counts in a quartet.
	size_t *columns_;				// 5|6|7|?
	Locus *site_;					// a vector to store the calls from reads.
	size_t *samples_;				// the number of samples (i.e. different individuals or populations) in the profile.
	std::vector <float_t> sample_gof_;		// the number of samples (i.e. different individuals or populations) in the profile.
public:
	const float_t getsample_property(const count_t &) const;
	const std::string getsample_name(const count_t &) const;
	char control;					//!< a variable that controls switches in the binary read/write mode.
	profile_header();				//!< does not initilize . . .
	profile_header(profile *);			//!< initilizes profile_header.
	~profile_header(void);
	void init(profile *);				//!< initilizes profile_header.

	const id0_t encodeid0(const std::string &);
	const id1_t encodeid1(const std::string &);
	const char encodeextraid(const char &, const size_t &);

	const std::string decodeid0(const id0_t&);
	const std::string decodeid1(const id1_t&);
	const std::string decodeextraid(const char &, const size_t &);

	int readheader(std::istream *);			//!< reads the header of a profile. All profiles from v 2.0 and later should have headers.
	int readtailer(std::istream *);			//!< reads the tailer of a profile. All profiles from v 2.0 and later should have headers.
	int writeheader(std::ostream *);		//!< write the . . .
	int writetailer(std::ostream *);		//!< write the . . .

	void set_delim_column(const char&);		//!< the delimiter which seperates columns
	void set_delim_quartet(const char&);		//!< the delimiter which seperates columns

	int setsamples(const count_t&);			//!< set the number of samples in the profile (only called in write mode).
	int setcolumns(const count_t&);			//!< set the number of columns for reading and writing.
	int setcolumn_name(const count_t&, const std::string&);
	const std::string getcolumn_name(const count_t&) const;
	profile_header & operator=(const profile_header&); //!< I don't think the copy makes sense. . . 
	
	void clear(void);
};

	
/*! \breif A table to store quartet information. 
 *	Eventually this will be a specialization to add conversion functions from mpileup, and the various .pro 
 *	format to a more map format that will be used to read/write information within mapgd commands. 
 *	See the design document for the map format in the docs directory. 
 */
class  profile{
//private variables should be initialized by reading the header...
private:

	/*these should all be controlled through the header_*/

	char delim_column;			//!< the delimiter which seperates columns
	char delim_quartet;			//!< the delimiter that seperates counts in a quartet
	size_t columns_;			//!< 5|6|7|more?

	size_t samples_;			//!< the number of samples (i.e. different individuals or populations) in the profile.
	id1_t size_;				//!< the number of lines in the sample. 0 if unkown.

	bool open_;				//!< indicates whether the profile opened succesfully
	bool sig_;				//!< were alleles thrown out if the allele only occurred in reads from one direction?
	bool read_;				//!< file is open for reading.
	bool write_;				//!< file is open for writing.
	bool binary_;				//!< binary mode flag. Incompatable with mpileup and noheader flags.
	bool mpileup_;				//!< file is an mpileup. Setting this should set the noheader flag.
	bool noheader_;				//!< file has no header

	/*done*/

	bool donothing_;			// a flag to indicate that nothing should be read for infile stream when read is called 

	static const count_t defaultorder[5];	// 01234 

	Locus site_;				//!< a structure that stores quartet information.
//	Locus site_buffer[BUFFER_SIZE];		//!< a structure that stores quartet information.

	int readm(Locus &);			//!< read file in mpileup mode.
	int readt(Locus &);			//!< read file in text mode.
	int readb(Locus &);			//!< read file in binary mode.

	int writet();				//!< write the quartet information in memory to file in text mode.
	int writeb();				//!< write the quartet information in memory to file in binary mode.
						// pro files do not contain enough information to construct mpileup files.

	int writet(Locus const &);		//!< write the quartet information passed to file in text mode.
	int writeb(Locus const &);		//!< write the quartet information passed to file in binary mode.

	profile_header header_;

	std::istream *in;			//!< all data is read from in.
	std::ostream *out;			//!< all data is writen is writen to out.
	std::fstream inFile;			//!< the file to read data from (if not stdin).
	std::ofstream outFile;			//!< the file to write data to (if not stdout).

	friend profile_header::profile_header(profile *);//!< Should set up pointers etc.
	friend void profile_header::init(profile *);	//!< Should set up pointers etc.
	void inline scan(const Locus &,const std::string &, quartet_t &); //?
public:
	profile();					//!< default constructor
	static const std::string names_;		//!< ACGTN
	const float_t getsample_property(const count_t &) const;

	profile* open(const char *, const char *);	//!< The function that opens a profile (if file).
	profile* open(const char *);			//!< The function that opens a profile (if stdin).
	bool is_open(void) const;			//!< Returns true if profile is open, false otherwise.

	/*basic io operation*/
	int copy(const profile&);			//!< Copys a line from profile
	int read();					//!< Reads a line from the instream. Returns 0 on success, EOF on EOF.
	int read(Locus &);				//!< Peaks at a line in the instream. Returns 0 on success, EOF on EOF.
	int write();					//!< Writes a line to the outstream. Returns 0 on success, EOF on EOF.
	int write(Locus const &);			//!< Writes a line to the outstream. Returns 0 on success, EOF on EOF.

							//!< Both of these functions should thow an error if the are used while streams are not open.
	int seek(const std::streampos);			//!< goes to the pos streampos of the stream. Returns 0 on success, EOF on if streampos is not in the stream. 
	int seek(std::string, std::string);		//!< goes to the line specified by the ID0 ID1 pair om the instream. Returns 0 on success, EOF if ID0 ID1 is not in the stream. 

	void close(void);				//!< close iostreams, writes tail, etc.

	void static merge(std::list <profile *>);	//!< take multiple profile and combine them into a single profile.

	/*functions dealing with the header*/

	int copyheader(const profile&);			//!< copys the header of a profile.
	int readheader();				//!< reads the header of a profile. All profiles from v 0.2.2 and later should have headers.
	int writeheader();				//!< Writes the header of a profile. All profiles from v 0.2.2 and later should have headers.
	
	int writetailer();				//!< Writes the header of a profile. All profiles from v 0.2.2 and later should have headers.

	void set_delim_column(const char&);		//!< The delimiter which seperates columns
	void set_delim_quartet(const char&);		//!< The delimiter which seperates columns

	void setsamples(count_t);			//!< Set the number of samples in the profile (only called in write mode).
	void setcolumns(count_t);			//!< Set the number of columns for reading and writing.
	size_t get_columns(void) const {return columns_;};	//!< Set the number of columns for reading and writing.

	int setcolumn_name(const count_t&, const std::string &);
	int setsample_name(const count_t&, const std::string &);

	const std::string getsample_name(const count_t &) const;

	/*a set of functions for converting the string information read in text mode into count_t and vice versa.*/

	const id0_t encodeid0(const std::string &);
	const id1_t encodeid1(const std::string &);
	const char encodeextraid(const char &, const count_t &);

	const std::string decodeid0(const id0_t &);
	const std::string decodeid1(const id1_t &);
	const std::string decodeextraid(const char &, const size_t &);

	const id0_t get_id0(void) const;
	const id1_t get_id1(void) const;
	const char get_extraid(const size_t &) const;

	void set_id0(const id0_t &);
	void set_id1(const id1_t &);
	void set_extraid(const char &, const size_t &);

	/*functions dealing with ?*/
	size_t size(void) const;			//!< Returns the number of samples/individuals in a profile.
	size_t length(void) const {return size_;};	//!< Returns the number of rows in a profile.

	/*functions dealing with the quartets*/

	void setbasecount(const count_t &, const count_t &, const count_t &);	//!< Sets sample x, base y to z.

	/*playing with the mask*/

	void mpileup(void) {mpileup_=true; noheader_=true;}	//!< File is an mpileup. Setting this should set the noheader flag.
	void noheader(void) {noheader_=true;}			//!< File has no header

	const id1_t getlinenumber(void) const;


	std::string getids(const Locus &);

	std::string getids(void);

	count_t maskedcount(void);			//!< Returns a count of the number of quartets masked at the current Locus. 
	void maskall(void);				//!< Mask all lines at the current Locus. 


	void header(void) {noheader_=false;};		//!< If header is true pro-file will expect a file to have a header...

	//locus

	Locus & get_locus(void) {return site_;};	//!< Returns the current Locus.
	void set_locus(const Locus &locus) {site_=locus;};	//!< Set the current Locus to locus.
	void sort(void);				//!< Left sort the sum of unmaksed quartets (i.e. record 0 is the larges).

	//quartet

	quartet_t get_quartet(const size_t &t) {return site_.get_quartet(t);};	//!< Returns the quartet of sample t at the current locus.
	std::vector <quartet_t>::iterator begin(void) {return site_.begin();};	//!< Returns an iterator to the quartet at the begining of the current locus.
	std::vector <quartet_t>::iterator end(void) {return site_.end();};	//!< Returns an iterator to the quartet just past the end of the current locus. */

};
	
	
#endif
