#ifndef LOCUS_H_
#define LOCUS_H_

#include <iostream>
#include <vector>
#include "../typedef.h"

#include "quartet.h"
#include "key.h"

/// A class to sort the alleles at [quartet_t]s.
/* The class locus handles sorting alleles at [quartet_t]s. Operations that use the relative count of [quartet_t]s etc.
 * should be members of this class.
 */ 
class locus {
private:	
	gt_t sorted_[5];								//!< an array to allow sorted access to quartets.
	std::vector <quartet_t> sample_;						//!< the five bases A/C/G/T/N;
public:
	size_t quartet_size(void) const {return sample_.size();};				//!< returns the number of quartet_t stored at locus.

	static constexpr gt_t defaultorder[5] = {0,1,2,3,4};				//!<???

	locus (const size_t &);
	locus ();

	/*BASIC OPERATORS*/
	locus & operator=(const locus&);		//!<
	locus & operator+=(const locus&);		//!<
	const locus operator+(const locus&) const;	//!<	

	gt_t get_index(const gt_t &a) const;		//!< returns the index of the a'th alleles in a sorted quartet_t.
	count_t get_coverage(const size_t &s) const;	//!< returns coverage of population/individual s.
	count_t get_coverage(void) const;		//!< returns total coverage.

	count_t get_count(const gt_t &a) const;		//!< returns the count of allele in all individuals in the population.
	count_t get_count(const size_t &s, const gt_t &a) const;//!< returns the count of individual/population s's a'th allele.

	void swap(const gt_t &lhs, const gt_t &rhs);	//!< exchange the alleles.
	void sort(const size_t &s);				//!< sort counts from most common to least common (based on poulation N).
	void sort(void);				//!< sort counts from most common to least common (among all non-masked sites).

	/**/
	void set_quartet(const quartet_t &q, const size_t &s);	//!< sets quartet_t of individual/population s to q (unsorted).

	const quartet_t & get_quartet(const size_t &s) const; 	//!< returns the s'th quartet at the locus.	
	quartet_t & get_quartet(const size_t &s);		//!< returns the s'th quartet at the locus.

	/*MASKING*/
	size_t maskedcount(void) const;			//!< returns the count of the number of individuals/populations that are masked.

	void maskall(void);					//!< mask all individuals/populations.
	void unmaskall(void);					//!< unmask all lines
	void mask(const std::vector <size_t> &v);		//!< mask all lines
	void unmask(const std::vector <size_t> &v);		//!< unmask all lines

	char get_name( const gt_t &a) const;			//!< get the nucleotide represented by the count at the a'th indexed quartet_t.
	char get_name_gt( const gt_t &a) const;			//!< get the nucleotide represented by the count at the a'th indexed quartet_t, return * if the count equals the N+1 indexed quartet.

	void resize(const size_t &);				//!< change the number of quartet_t s at the locus.

	void mask_low_cov(const count_t &dp);						//!< mask site with coverage strictly less than dp;

	std::vector <quartet_t>::iterator begin(void) {return sample_.begin();};	//!< return an iterator to the quartet_t s stored at this locus.
	std::vector <quartet_t>::iterator end(void) {return sample_.end();};		//!< return an iterator to the quartet_t s stored at this locus.

	std::vector <quartet_t>::const_iterator begin(void) const {return sample_.begin();};	//!< return an iterator to the quartet_t s stored at this locus.
	std::vector <quartet_t>::const_iterator end(void) const {return sample_.end();};		//!< return an iterator to the quartet_t s stored at this locus.

	//ostream& operator<<(ostream& output, const Point& p);
};

size_t size_bin(const locus &l) {return l.quartet_size()*sizeof(quartet_t)+sizeof(gt_t)*5;};	//!< all class of type data need to decleare a size function that returns the total size of ...
size_t size_tex(const locus &l) {return l.quartet_size()*4;};	//!< all class of type data need to decleare a size function that returns the total size of ...
#endif
