#ifndef QUARTET_H_
#define QUARTET_H_

#include "typedef.h"
#include <sstream>
#include <cstring> 
#include <iostream>
#include <iomanip> 

/// A class that stores quartet information.
/*! 
 *
 */
typedef struct quartet {
	count_t base[5];	//!< The count of the number of occurnaces of bases. What nucleotieds this cout represents is stored at the Locus. 
	bool masked;		//!< A flag that indicates whether or not functions should use information from this quartet_t. Functions may ignore this flag if such behavior is desired. 
	char delim;
	
	/*! The default constructor unmaskeds the site and zeros the base counts.
	 */
	quartet (){		
		masked=false;	
		memset(base, 0, 5*sizeof(count_t) );
		delim='/';
	}
	
	/*! This constructor explicitly sets the counts of A, C, G, T, and N and unmasked the site.
	 */
	quartet(const count_t &A, const count_t &C, const count_t &G, const count_t &T, const count_t &N) {
		base[0]=A;
		base[1]=C;
		base[2]=G;
		base[3]=T;
		base[4]=N;
		masked=false;
	}

	/*! This operator is not implemented correctly. TODO fix it.
	 */
	quartet& operator+=(const quartet& x) {
		memcpy(base, x.base, 5*sizeof(count_t) );
		return *this;
	}

	/*! This operator adds the basecounts of the two quartet s together and returns them. masked is ignored.
	 */
	inline quartet operator+(const quartet& x) const {
		return quartet(base[0]+x.base[0], base[1]+x.base[1], base[2]+x.base[2], base[3]+x.base[3], base[4]+x.base[4]);
	}

	/*! This operator copies the array count_t to base. Mask and delimiter is left untouched.
	 */
	quartet& operator=(const count_t& x) {
		memset (base, x, 5*sizeof(count_t) );
		return *this;
	}
} quartet_t;

/*! \breif returns the quartet */
count_t count(const quartet_t&);

/*! \breif Sets masked to true.*/
void mask(quartet_t&);

/*! \breif returns the */
void unmask(quartet_t&);

//writet(stream, quartet_t&);
//writet(stream, quartet_t&);

std::ostream& operator<< (std::ostream& out, const quartet& q);
std::istream& operator>> (std::istream& in, quartet& q);
#endif
