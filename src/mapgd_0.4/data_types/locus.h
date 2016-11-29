#ifndef LOCUS_H_
#define LOCUS_H_

#include <iostream>
#include <vector>

#include "quartet.h"
#include "base.h"
#include "stream-tools.h"
#include "data.h"
#include "file-index.h"

/* This is the record stored in a pro file*/
class Locus : virtual public Indexed_data{
private:
	//NOT READ READ/WRITE!!
	std::vector <std::string> sample_names_;		//!< names of the samples sequenced.

	//READ READ/WRITE!!
	using Indexed_data::abs_pos_;

	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Locus(Columns);
	}
	void write (std::ostream& out) const;
	void read (std::istream& in);
public:
	void write_binary (std::ostream& out) const;
	void read_binary (std::istream& in);
	/* BEGIN DATA BLOCK */
	/* these need to be changed to private */
	gt_t sorted_[5];				//!< an array to allow sorted access to quartets.
	std::vector <quartet_t> sample;			//!< The five bases A/C/G/T/N;
	Base ref;
	/* END DATA BLOCK */

	static const gt_t default_order[5];
//public:

	Locus (const count_t &);
	Locus ();
	Locus (const std::vector<std::string> &); 

	/*BASIC OPERATORS*/
	Locus & operator=(const Locus&);	
	Locus & operator+=(const Locus&);	
	const Locus operator+(const Locus&) const;	
	/** }*/

	const count_t getindex(count_t) const;		//!< Returns the index of the alleles in order a sorted order.

	count_t getcoverage(count_t) const;		//!< Returns coverage of population/individual N.
	count_t getcoverage(void) const;		//!< Returns total coverage.

	count_t getcount(count_t) const;		//!< Returns the count of allele in all individuals in the population.
	count_t getcount(count_t, count_t) const;	//!< Returns the count of individual a's b'th allele.

	void swap(count_t, count_t);			//!< Exchange the alleles.
	void sort(count_t);				//!< Sort counts from most common to least common (based on population N).
	void sort(void);				//!< Sort counts from most common to least common (among all non-masked sites).

	/* \defgroup QUARTET Quartet
	 * The set of functions dealing with quartets
	 * @{
	 */
	void set_quartet(const quartet_t &, const count_t &);	//!< Sets the quartet_t array (unsorted).
	const quartet_t & get_quartet(const count_t &) const; 	//!< Returns the N'th quartet at the Locus.	
	quartet_t & get_quartet(const count_t &);		//!< Returns the N'th quartet at the Locus.
	 /** @} */

	/* \defgroup MASKING Masking 
	 * Functions which set or unset the masked flag, quartets which are 
	 * masked are ignored by calculations, and are printed as 0's. 
	 * @{
	 */
	count_t maskedcount(void) const;		//!< Returns the count of the number of individuals that are masked.
	void maskall(void);				//!< Mask all lines
	void unmaskall(void);				//!< Unmask all lines
	void mask(const std::vector <size_t> &);	//!< Mask all lines
	void unmask(const std::vector <size_t> &);	//!< Unmask all lines
	void mask_low_cov(const count_t &dp);	//!< Mask site with coverage stricktly lt dp;
	/** @} */

	/* \defgroup DEPRICATED depricated 
	 * Functions that were once used, but are no longer usful.
	 * @{
	 */
		char getname( const count_t &) const;			//!< Get the nucleotide represented by the count at the N'th indexed quartet_t.
		char getname_gt( const count_t &) const;		//!< Get the nucleotide represented by the count at the N'th indexed quartet_t, return * if the count equals the N+1 indexed quartet.
		count_t get_extraid(const size_t &) const;		//!< Get the extraid of the Locus. Used to represent the reference call.
		void set_extraid(const count_t &, const size_t &);	//!< Set the extraid of the Locus. Just used to represent the reference call. 
	/* this will definetly by dropped */
	/** @} i
	  */


	using Indexed_data::get_abs_pos;		//!< Get the absolute position of the Locus.
	using Indexed_data::set_abs_pos;		//!< Set the absolute position of the Locus.

	void resize(const size_t &);			//!< Change the number of quartet_t s at the Locus.

	/* \defgroup iterators Iterators
	 * Locus can return the iterators foir the quartets.
	 * @{
	 */
	inline std::vector <quartet_t>::iterator begin(void) {return sample.begin();};		//!< Return an iterator to the quartet_t s stored at this Locus.

	//! Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::iterator 
	end(void) {
		return sample.end();
	};

	//! Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::const_iterator 
	cbegin(void) const {
		return sample.cbegin();
	};		

	//! Return an iterator to the quartet_t s stored at this Locus.
	inline std::vector <quartet_t>::const_iterator 
	cend(void) const {
		return sample.cend();
	};		
	/** @}*/

	//! names of the samples sequenced
	inline std::vector <std::string> 
	get_sample_names(void) const 
	{
		return sample_names_;
	};

	//! names of the samples sequenced.
	inline void 
	set_sample_names(const std::vector <std::string>& sample_names) 
	{
		sample_names_=sample_names;
		if (sample_names_.size()!=sample.size() ) sample.assign(sample_names_.size(), quartet() );
	};		

	/* \defgroup LOCUS_DATA Inherited members
	 * These members are inherited from Data
	 * @{
	 */
	std::string header(void) const;
	size_t size(void) const;

	static const std::string file_name;
	static const std::string table_name;
	static const bool binary;

	const std::string get_file_name(void) const;
	const std::string get_table_name(void) const;
	const bool get_binary(void) const;

	const std::string sql_header(void) const;				
	const std::string sql_column_names(void) const;				
	const std::string sql_values(void) const;				

	void sql_read(std::istream &) override;
	/** @}*/
	
	friend std::istream& mpileup (std::istream& in, Locus& x, const int &, const int &);
	friend void scan(const Locus & site, const std::string &str, quartet_t &q);
};


#endif
