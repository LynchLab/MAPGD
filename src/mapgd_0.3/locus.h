#ifndef LOCUS_H_
#define LOCUS_H_

#include <iostream>
#include <vector>
#include "quartet.h"

class Locus {
private:
public:
	/* these need to be changed to private */
	count_t sorted_[5];				// an array to allow sorted access to quartets.
	
	std::vector <quartet_t> sample;			//The five bases A/C/G/T/N;
	std::vector <count_t> extraid;			//extra ids associated with the quartet. (ref base identiy?).

//TODO count_t to be replaced.
	id0_t id0;
	id1_t id1;					//The ids associated with the quartet.
//public:

	Locus (count_t);
	Locus ();

	/*BASIC OPERATORS*/
	Locus & operator=(const Locus&);	
	Locus & operator+=(const Locus&);	
	const Locus operator+(const Locus&) const;	



	const count_t getindex(count_t) const;		//returns the index of the alleles in order a sorted order

	count_t getcoverage(count_t) const;		//returns coverage of population/individual N
	count_t getcoverage(void) const;		//returns total coverage

	count_t getcount(count_t) const;		//returns the count of allele in all indivuals in the population.
	count_t getcount(count_t, count_t) const;	//returns the count of individual a's b'th allele.

	void swap(count_t, count_t);			//exchage the alleles
	void sort(count_t);				//sort reads from most common to least common (based on poulation N).
	void sort(void);				//sort reads from most common to least common (amoung all non-masked sites).

	/**/
	void set_quartet(const quartet_t &, const count_t &);	//sets the quartet array (unsorted)

	const quartet_t & get_quartet(const count_t &) const;	
	quartet_t & get_quartet(const count_t &);	

	/*MASKING*/
	count_t maskedcount(void) const;		//returns the count of the number of individuals that are masked.

	void maskall(void);				//mask all lines
	void unmaskall(void);				//unmask all lines

	char getname( const count_t &) const;
	char getname_gt( const count_t &) const;

	id0_t get_id0(void) {return id0;};
	id1_t get_id1(void) {return id1;};

	void set_id0(const id0_t &tid0) {id0=tid0;};
	void set_id1(const id1_t &tid1) {id1=tid1;};
	void set_extraid(const count_t &, const size_t &);

	count_t get_extraid(const size_t &) const;

	void resize(const size_t &);

	std::vector <quartet_t>::iterator begin(void) {return sample.begin();};	
	std::vector <quartet_t>::iterator end(void) {return sample.end();};
};

#endif
