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
#include "typedef.h"
#include "streamtools.h"

#define NONE		0
#define BADHEADER	1
#define UNEXPECTED	2
#define VER		2.0

#define SKIP		1

typedef char name_t;

typedef struct quartet {
	count_t base[5];
	bool masked;
	
	quartet (){
		masked=false;
		memset (base,0, 4*sizeof(count_t) );
	}
	
	quartet(count_t A, count_t C, count_t G, count_t T, count_t N){
		base[0]=A;
		base[1]=C;
		base[2]=G;
		base[3]=T;
		base[4]=N;
	}

	quartet& operator+=(const quartet& rhs)
	{
		base[0]+=rhs.base[0];
		base[1]+=rhs.base[1];
		base[2]+=rhs.base[2];
		base[3]+=rhs.base[3];
		base[4]+=rhs.base[4];
		return *this;
	}
	inline quartet operator+(const quartet& rhs) const
	{
		return quartet(base[0]+rhs.base[0], base[1]+rhs.base[1], base[2]+rhs.base[2], base[3]+rhs.base[3], base[4]+rhs.base[4]);
	}

} quartet_t;

count_t major(const quartet_t);
count_t minor(const quartet_t);
count_t count(const quartet_t);

class site_t {
private:
public:
	site_t (count_t);
	site_t ();
	std::vector <quartet_t> sample;			//The five bases A/C/G/T/N;
	std::string id0, id1;				//The ids associated with the quartet.
	std::vector <std::string> extra_ids;		//extra ids associated with the quartet. (ref base identiy?).
	site_t & operator=(const site_t&);	
};

class  profile{

private:
	bool open_;					// indicates whether the profile opened succesfully
	char delim_column='\t';				// the delimiter which seperates columns
	char delim_quartet='/';				// the delimiter that seperates counts in a quartet

	unsigned int columns_;				// 5|6|7|more?
	unsigned int samples_;				// the number of samples (i.e. different individuals or populations) in the profile.
	double sig_;					// were alleles thrown out if the allele only occurred in reads from one direction?
	site_t site_;					// a vector to store the calls from reads
	count_t sorted_[5];				// an array to sort reads.

	bool read_;					// file mode flag;
	bool write_;					// file mode flag;

	static const std::string names_;		// the names of the bases : ACGTN.

	bool donothing_;				// a flag to indicate that nothing should be read for infile stream when read is called 
	std::vector <std::string> column_names;		// the names of all the columns in the profile.
							// (this will be depricated when headers are finished).

	std::istream *in;				// all data is read from in.
	std::ostream *out;				// all data is writen is writen to out.
	std::fstream inFile;				// the file to read data from (if not stdin).
	std::ofstream outFile;				// the file to write data to (if not stdout).

public:
	profile();					//default constructor

	profile* open(const char*, const char);		//The function that opens a profile (if file).
	profile* open(const char);			//The function that opens a profile (if stdin).
	bool is_open(void) const;			//returns true if profile is open, false otherwise.

	int copy(const profile&);			//copys a line from profile
	int copyheader(const profile&);			//copys the header of a profile.

	int read();					//reads a line from the instream. Returns 0 on success, EOF on EOF.
	int readheader();				//reads the header of a profile. All profiles from v 2.0 and later should have headers.
							//Returns 0 on success, EOF on EOF (should never happen), and a failure code on faliure?.
	int read(int);					//peaks at a line in the instream. Returns 0 on success, EOF on EOF.
	int write();					//writes a line to the outstream. Returns 0 on success, EOF on EOF.
	int write(site_t const &);			//writes a line to the outstream. Returns 0 on success, EOF on EOF.
	int writeheader();				//reads the header of a profile. All profiles from v 2.0 and later should have headers.

	int seek(const std::streampos);//TODO: Implement.		//goes to line N of the instream. Returns 0 on success, EOF on EOF. 
	int seek(std::string, std::string);//TODO: Implement.		//goes to line N of the instream. Returns 0 on success, EOF on EOF. 

	void sort(count_t);				//sort reads from most common to least common (based on poulation N).
	void sort(void);				//sort reads from most common to least common (amoung all non-masked sites).

	void set_delim_column(const char&);		// the delimiter which seperates columns
	void set_delim_quartet(const char&);		// the delimiter which seperates columns

	count_t getcount(count_t) const;		//returns the population count.
	const count_t *getquartet(count_t) const;	//returns the quartet array (unsorted)

	std::vector <quartet_t>::const_iterator begin(void) const;	//returns a pointer to the begining (unsorted)
	std::vector <quartet_t>::const_iterator end(void) const;	//returns a pointer to the end (unsorted)

	std::vector <quartet_t>::iterator begin(void);	//returns a pointer to the begining (unsorted)
	std::vector <quartet_t>::iterator end(void);	//returns a pointer to the end (unsorted)

	count_t getcount(count_t, count_t) const;	//returns the count of individuals a's b'th allele.

	count_t maskedcount(void) const;		//returns the count of the number of individuals that are masked.

	name_t getname(count_t) const;			//returns the name [*i.e. ACG or T] of the sorted alleles.
	name_t getname_gt(count_t) const;			//returns the name [*i.e. ACG or T] of the sorted alleles.
	const count_t getindex(count_t) const;		//returns the index of the alleles in order a sorted order
	std::string getids(void) const; 		//
	count_t getcoverage(count_t) const;		//returns coverage of population/individual N
	count_t getcoverage(void) const;		//returns total coverage
	count_t size(void) const;			//number of populations/individuals
	void maskall(void);				//mask all lines
	void unmask(count_t);				//unmask line N
	void unmask(quartet *);				//mask line N
	void mask(quartet *);				//mask line N
	void mask(count_t);				//mask line N
	int setsamples(count_t);			//set the number of samples in the profile (only called in write mode).
	int setcolumns(count_t);			//set the number of columns for reading and writing.
	int setcolumn_name(count_t, std::string str);

	void close(void);				//close iostreams.
	void swap(count_t, count_t);
	void static merge(std::list <profile *>);	//take multiple profile and combine them into a single profile.
};
	
#endif
