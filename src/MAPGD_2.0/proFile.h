#ifndef PROFILE_HPP_
#define PROFILE_HPP_	

#include <cstdio>
#include <cstring>
#include <iostream>
#include <climits>

#include <list>
#include <vector>
#include <algorithm>
#include "streamtools.h"

typedef unsigned int count_t;
typedef char name_t;
typedef struct quartet {
	count_t base[5];
} quartet_t;

typedef struct site {
	std::vector <quartet_t> sample;			//The five bases A/C/G/T/N;
	std::string id0, id1;				//The ids associated with the quartet.
} site_t;

class  profile{

private:
	bool open_;					// indicates whether the profile opened succesfully
	char delim='\t';				// the feild delimiter
	unsigned int columns_;				// 5|6?
	unsigned int samples_;				// 5|6?
	double sig_;					// were alleles thrown out if read from only one direction?
	site_t site_;					// a vector to store the calls from reads
	count_t sorted_[5];				// a vector to store the calls from reads
	static const std::string names_;		// the names of the bases ACGTN.
	bool *masked_;					// an array to determin which populations to use for group calculations.	

	std::istream *in;				// all data is read from in.
	std::ostream *out;				// all data is writen is writen to out.
	std::ifstream inFile;				// the file to read data from (if not stdin).
	std::ofstream outFile;				// the file to write data to (if not stdout).

public:
	profile();					//default constructor

	profile* open(const char*, const char);		//The function that opens a profile.
	profile* open(const char);			//The function that opens a profile.
	bool is_open(void);				//returns true if profile is open, false otherwise.

	int read(void);					//reads a line from the instream. Returns 0 on success, EOF on EOF.
	int write(void);				//reads a line from the instream. Returns 0 on success, EOF on EOF.
	int seek(void);					//reads a line from the instream. Returns 0 on success, EOF on EOF.

	void sort(count_t);				//sort reads from most common to least common.
	void sort(void);				//sort reads from most common to least common.

	count_t getcount(count_t);			//
	count_t getcount(count_t, count_t);		//
	name_t getname(count_t);			//
	std::string getids(void);			//
	count_t getcoverage(count_t);			//
	count_t getcoverage(void);			//
	count_t samplesize(void);			//
	void maskall(void);				//
	void unmask(count_t);				//
	void mask(count_t);				//

	void close(void);				//close instream.
	void static merge(std::list <profile *>);	//take multiple profile and combine them into a single profile.
};
	

#endif
