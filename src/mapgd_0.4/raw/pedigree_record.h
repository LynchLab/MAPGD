#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_

#include "typedef.h"
#include <iostream>

class Pedigree_record {
private :
	Pedigree *pedigree_;
public :
	std::string name;
	uint8_t sex;
	size_t family;						
	std::vector <Pedigree_record *> parent, offspring;
	Pedigree_record( const std::vector <Pedigree_record *> &, const std::vector <Pedigree_record *> &, Pedigree *);	//!< constructor.
	Pedigree_record( const std::vector <Pedigree_record *> &, Pedigree *, const bool &);				//!< constructor.
	Pedigree_record( Pedigree *);											//!< constructor.
	Pedigree_record();												//!< constructor.
	Pedigree_record & operator= (const Pedigree_record&);
	friend std::ostream& operator << (std::ostream& out, const Pedigree_record& x);
	friend std::istream& operator >> (std::istream& in, Pedigree_record& x);
};

#endif
