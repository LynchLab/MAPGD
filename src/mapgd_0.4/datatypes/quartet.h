#ifndef QUARTET_H_
#define QUARTET_H_

#include "../typedef.h"
#include <cstring> 

/// A structure that stores quartet information.
typedef struct quartet {
	count_t base[5];	//!< The count of the number of occurnaces of bases. 
				//What nucleotieds this cout represents is stored at the locus. 
	bool masked;		//!< A flag that indicates whether or not functions should use information from this quartet_t.
				// Functions may ignore this flag if such behavior is desired. 
	
	/// The default constructor unmasks the site and zeros the base counts.
	quartet (){		
		masked=false;	
		memset(base, 0, 5*sizeof(count_t) );
	}
	
	/// This constructor explicitly sets the counts of A, C, G, T, and N and unmasked the site.
	quartet(const count_t &A, const count_t &C, const count_t &G, const count_t &T, const count_t &N) {
		base[0]=A;
		base[1]=C;
		base[2]=G;
		base[3]=T;
		base[4]=N;
		masked=false;
	}

	/// This operator is not implemented correctly. TODO fix it.
	quartet& operator+=(const quartet& x) {
		memcpy(base, x.base, 5*sizeof(count_t) );
		return *this;
	}

	/// This operator adds the basecounts of the two quartet s together and returns them. masked is ignored.
	inline quartet operator+(const quartet& x) const {
		return quartet(base[0]+x.base[0], base[1]+x.base[1], base[2]+x.base[2], base[3]+x.base[3], base[4]+x.base[4]);
	}

	/// This operator copies the array count_t to base. mask is left untouched.
	quartet& operator=(const count_t& x) {
		memset (base, x, 5*sizeof(count_t) );
		return *this;
	}
} __attribute__((packed)) quartet_t;

/// returns the total depth of coverge at the quartet
count_t count(const quartet_t&);

/// sets masked to true
void mask(quartet_t&);

/// sets masked to false
void unmask(quartet_t&);

#endif
