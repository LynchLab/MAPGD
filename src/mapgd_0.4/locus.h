#ifndef LOCUS_H_
#define LOCUS_H_

#include <iostream>
#include <vector>
#include "quartet.h"
#include "base.h"
#include "stream-tools.h"
//#include "pro-file.h"

/* This is the record stored in a pro file*/

class Locus {
private:
	std::vector <std::string> sample_names_;		//!< names of the samples sequenced.
public:
	/* these need to be changed to private */
	gt_t sorted_[5];				//!< an array to allow sorted access to quartets.
	
	std::vector <quartet_t> sample;			//!< The five bases A/C/G/T/N;
	std::vector <count_t> extraid;			//!< extra ids associated with the quartet. (ref base identiy?).

	id0_t id0;
	id1_t id1;					//!< The ids associated with the quartet.
	Base ref;
	static const gt_t default_order[5];
//public:

	Locus (const count_t &);
	Locus ();
	Locus (const std::vector<std::string> &); 

	/*BASIC OPERATORS*/
	Locus & operator=(const Locus&);	
	Locus & operator+=(const Locus&);	
	const Locus operator+(const Locus&) const;	

	const count_t getindex(count_t) const;		//!< Returns the index of the alleles in order a sorted order.

	count_t getcoverage(count_t) const;		//!< Returns coverage of population/individual N.
	count_t getcoverage(void) const;		//!< Returns total coverage.

	count_t getcount(count_t) const;		//!< Returns the count of allele in all indivuals in the population.
	count_t getcount(count_t, count_t) const;	//!< Returns the count of individual a's b'th allele.

	void swap(count_t, count_t);			//!< Exchage the alleles.
	void sort(count_t);				//!< Sort counts from most common to least common (based on poulation N).
	void sort(void);				//!< Sort counts from most common to least common (amoung all non-masked sites).

	/**/
	void set_quartet(const quartet_t &, const count_t &);	//!< Sets the quartet_t array (unsorted).

	const quartet_t & get_quartet(const count_t &) const; 	//!< Returns the N'th quartet at the Locus.	
	quartet_t & get_quartet(const count_t &);		//!< Returns the N'th quartet at the Locus.

	/*MASKING*/
	count_t maskedcount(void) const;		//!< Returns the count of the number of individuals that are masked.

	void maskall(void);				//!< Mask all lines
	void unmaskall(void);				//!< Unmask all lines
	void mask(const std::vector <size_t> &);	//!< Mask all lines
	void unmask(const std::vector <size_t> &);	//!< Unmask all lines


	/* these will probably be dropped once I get ride of the byzantian pro-file*/
	char getname( const count_t &) const;		//!< Get the nucleotide represented by the count at the N'th indexed quartet_t.
	char getname_gt( const count_t &) const;	//!< Get the nucleotide represented by the count at the N'th indexed quartet_t, return * if the count equals the N+1 indexed quartet.

	/* dumb dumb dumb dumb */
	id0_t get_id0(void) const {return id0;};	//!< Get the id0 of the Locus.
	id1_t get_id1(void) const {return id1;};	//!< Get the id1 of the Locus.

	void set_id0(const id0_t &tid0) {id0=tid0;};	//!< Set the id0 of the Locus.
	void set_id1(const id1_t &tid1) {id1=tid1;};	//!< Set the id1 of the Locus.

	/* this will definetly by dropped */
	void set_extraid(const count_t &, const size_t &);	//!< Set the extraid of the Locus. Just used to represent the reference call. 

	count_t get_extraid(const size_t &) const;	//!< Get the extraid of the Locus. Used to represent the reference call.

	void resize(const size_t &);			//!< Change the number of quartet_t s at the Locus.

	/* inline functions */
	inline std::vector <quartet_t>::iterator begin(void) {return sample.begin();};		//!< Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::iterator end(void) {return sample.end();};		//!< Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::const_iterator cbegin(void) const {return sample.cbegin();};		//!< Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::const_iterator cend(void) const {return sample.cend();};		//!< Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <std::string> get_sample_names(void) const {return sample_names_;};		//!< names of the samples sequenced.
	inline void set_sample_names(const std::vector <std::string>& sample_names) {
		sample_names_=sample_names;
		if (sample_names_.size()!=sample.size() ) sample.assign(sample_names_.size(), quartet() );
	};		//!< names of the samples sequenced.

	void mask_low_cov(const count_t &dp);	//!< Mask site with coverage stricktly lt dp;

	/* group obligate functions for i/o */
	std::string header(void) const;
	size_t size(void) const;

	static const std::string file_name;
	static const std::string table_name;

	friend  std::ostream& operator<< (std::ostream& out, const Locus& x);
	friend std::istream& operator>> (std::istream& in, Locus& x);
	friend std::istream& mpileup (std::istream& in, Locus& x);
	friend void scan(const Locus & site, const std::string &str, quartet_t &q);
};


#endif
