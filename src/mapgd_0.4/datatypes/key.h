#ifndef _KEY_H_
#define _KEY_H_

#include <typeinfo>
#include <cstring>
#include <string>
#include <iostream>
#include <vector>
#include <map>

#include "../typedef.h"
#include "../stream_tools.h"
//#include "formats.h"

/// a key is a class that allows us to read/write data to rows. 
class key <template type>{
protected:
	uint8_t offset_;	//!< the offset of from the begining of a row
	type * value_;		//!< a pointer to the data this key stores
	size_t size_;		//!< The size of this object. 
	uint8_t instance_;	//!< if more than one instance of this key is stored in a row, this needs to 
				//!< be incrimented.

	key <char *> keyname_;	//!< the name of the key. See key_defs.txt. i
	key <keyid_t> keyid_;	//!< the number of the key. See key_defs.txt.
	ket <char *> keydesc_;	//!< the description of the key. 

	size_t (*size_of_f) (void *);
	void (*to_text_f) (uint8_t *, void*);
	void (*from_tx_f) (void *, uint8_t*);
public:
	static constexpr uint8_t nokey=-1;	//!< used to represent an uninstated key. 
	void set_offset(const uint8_t &new_offset);
	size_t size (void *);

	key(const char *init);		//!< returns a key ...

	key(const uint8_t &id);			//!< returns a key ... 

	key(void);				//!<

	const size_t& size(void) const {return size_;}			//!< returns size of the data referenced by this key in bytes.
	const uint8_t& offset(void) const {return offset_;}		//!< returns the number the offset (in bytes) from the beging of a row. 
	type & value (void);					//!< returns a pointer to the data refernecd by this key.

	const char * name(void) const;
	const uint8_t id(void) const;
	const char * description(void) const;

	key& operator=(const key &);
};


#endif
