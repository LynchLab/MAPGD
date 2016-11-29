#ifndef _REGION_H_
#define _REGION_H_

#include "typedef.h"
#include "stream-tools.h"
#include <iostream>
#include <vector>
#include "file-index.h"

class Region {
private :
public :
	id1_t abs_start, start, abs_stop, stop;
	id0_t id0;
	std::string scf_name;
	Region(const std::string &, const id1_t &, const id1_t &);	//!< constructor.
	Region(const id0_t &, const id1_t &, const id1_t &);		//!< constructor.
	Region(const id1_t &, const id1_t &);				//!< constructor.
	Region();	//!< constructor.
	friend std::ostream& operator << (std::ostream& out, const Region& x);
	friend std::istream& operator >> (std::istream& in, Region& x);
	void set(const File_index &);
};

bool isregion(const char *);
Region ator(const char *);

#endif
