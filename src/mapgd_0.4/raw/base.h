/* synonym for population? */

#ifndef BASE_H_
#define BASE_H_

#include <iostream>
#include "typedef.h"
#include <cfloat>
#include <iomanip>

/// A class that handels conversions between the human readable bases 
/// (A,C,G, and T) and their numberical representation.
/** 
 */
class Base { 
private:
	bool masked;
public:
	Base();
	Base(const char &);
	Base(const gt_t &);

	gt_t base;

	friend std::ostream& operator << (std::ostream&, const Base&);	//!< use the << operator to write Allele.
	friend std::istream& operator >> (std::istream&, Base&);		//!< use the >> operator to read Allele.

	static char btoc(const gt_t &);
	static gt_t ctob(const char &);
};

#endif
