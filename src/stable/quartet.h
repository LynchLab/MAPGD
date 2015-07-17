#ifndef QUARTET_H_
#define QUARTET_H_

#include "typedef.h"
#include <cstring> 

typedef struct quartet {
	count_t base[5];
	bool masked;
	
	quartet (){
		masked=false;
		memset (base,0, 5*sizeof(count_t) );
	}
	
	quartet(const count_t &A, const count_t &C, const count_t &G, const count_t &T, const count_t &N) {
		base[0]=A;
		base[1]=C;
		base[2]=G;
		base[3]=T;
		base[4]=N;
	}

	quartet& operator+=(const quartet& x) {
		memcpy(base, x.base, 5*sizeof(count_t) );
		return *this;
	}
	inline quartet operator+(const quartet& x) const {
		return quartet(base[0]+x.base[0], base[1]+x.base[1], base[2]+x.base[2], base[3]+x.base[3], base[4]+x.base[4]);
	}
	quartet& operator=(const count_t& x) {
		memset (base, x, 5*sizeof(count_t) );
		masked=false;
		return *this;
	}

	quartet& operator=(const quartet& rhs) {
		memcpy (base, rhs.base, 5*sizeof(count_t) );
		masked=rhs.masked;
		return *this;
	}
} quartet_t;

/* \breif returns the */
count_t major(const quartet_t&);
count_t minor(const quartet_t&);
count_t error1(const quartet_t&);
count_t error2(const quartet_t&);
count_t count(const quartet_t&);
#endif
