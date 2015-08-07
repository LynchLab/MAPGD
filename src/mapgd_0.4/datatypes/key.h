#ifndef _KEY_H_
#define _KEY_H_

#include <typeinfo>
#include <cstring>
#include <string>
#include <iostream>
#include "../typedef.h"
#include "basic_types.h"

/// a key is a class that allows for us to read/write data to rows. 
class key {
private:
	uint8_t row_offset_;		//!< the offset for storing this object in a row.
	void * value_;			//!<  
	size_t size_;			//!< The size of this object. 

	KEYNAME keyname_;		//!< the name of the key. See key_defs.txt. 
	KEYNUM keynum_;			//!< the number of the key. See key_defs.txt.
	KEYDESC keydesc_;		//!< the description of the key. See key_defs.txt.

public:
	void set_offset(const uint8_t &new_offset){offset_=new_offset;}; //!< should be made privateish.

	key(const std::string &name, 
		const uint8_t &num, 
		const std::string &desc, 
		void * data);			//!< returns a key ...
	
	static constexpr uint8_t nokey=-1;	//!< used to represent an uninstated key. 

	 /* The objects pointed to by keys can be containers, so sizeof() will not return the correct size of the object.
	  *  However, objects pointed to by keys must be allocated in single contigious block of memory of fixed size,
	  *  since I use memcpy to move stuff around. */

	uint8_t instance;			//!< if more than one instance of this key is stored in a row, this needs to 
						// be incrimented.

	size_t size(void) const {return size_;}			//!< returns size of the data referenced by this key in bytes.
	uint8_t offset(void) const {return offset_;}		//!< returns the number the offset (in bytes) from the beging of a row. 
	virtual void * value (void)=0;

	std::string get_name(void) const { return keyname_.value();};
	uint8_t get_num(void) const { return keynum_.value();};
	std::string get_desc(void) const { return keydesc_.value();};
};


namespace keyid{
	constexpr uint8_t genprob=16;
	constexpr uint8_t freq=17;
	constexpr uint8_t error=18;
	constexpr uint8_t loglike=19;
	constexpr uint8_t smpname=20;
	constexpr uint8_t rowid=21;
	constexpr uint8_t cov=22;
	constexpr uint8_t smpsize=23;
	constexpr uint8_t somatic=24;
	constexpr uint8_t valid=25;
	constexpr uint8_t strbias=26;
	constexpr uint8_t sb=26;
	constexpr uint8_t mapq0=27;
	constexpr uint8_t hapmap2=28;
	constexpr uint8_t end=29;
	constexpr uint8_t ref=30;
	constexpr uint8_t anctype=31;
	constexpr uint8_t ncut=32;
	constexpr uint8_t gof=33;
//	constexpr uint8_t default=34;
	constexpr uint8_t polyll=35;
	constexpr uint8_t hwell=36;
	constexpr uint8_t maxll=37;
	constexpr uint8_t efchrom=38;
}

#endif
