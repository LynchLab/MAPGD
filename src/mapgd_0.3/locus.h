#ifndef LOCUS_H_
#define LOCUS_H_

#include "quartet.h"

class site_t {
private:
public:
	count_t sorted_[5];				// an array to allow sorted access to quartets.
	site_t (count_t);
	site_t ();
	~site_t ();
	std::vector <quartet_t> sample;			//The five bases A/C/G/T/N;

	count_t samples_;

	count_t id0;
	uint64_t id1;				//The ids associated with the quartet.
	std::vector <count_t> extraid;		//extra ids associated with the quartet. (ref base identiy?).

	site_t & operator=(const site_t&);	

	const count_t getindex(count_t) const;		//returns the index of the alleles in order a sorted order
	count_t getcoverage(count_t) const;		//returns coverage of population/individual N
	count_t getcoverage(void) const;		//returns total coverage
	count_t getcount(count_t) const;		//returns the population count.
	count_t getcount(count_t, count_t) const;	//returns the count of individuals a's b'th allele.

	char getname(count_t) const;			//returns the name [*i.e. ACG or T] of the sorted alleles.
	char getname_gt(count_t) const;		//returns the name [*i.e. ACG or T] of the sorted alleles.

	void swap(count_t, count_t);			//exchage the alleles
	void sort(count_t);				//sort reads from most common to least common (based on poulation N).
	void sort(void);				//sort reads from most common to least common (amoung all non-masked sites).
	
	count_t maskedcount(void) const;		//returns the count of the number of individuals that are masked.
	const count_t *getquartet(count_t) const;	//returns the quartet array (unsorted)

	void maskall(void);				//mask all lines
	void unmaskall(void);				//mask all lines
	void unmask(count_t);				//unmask line N
	void unmask(quartet *);				//mask line N
	void mask(quartet *);				//mask line N
	void mask(count_t);				//mask line N
};

#endif
