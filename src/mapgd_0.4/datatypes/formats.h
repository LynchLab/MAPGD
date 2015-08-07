#ifndef _FORMAT_H_
#define _FORMAT_H_

#include <string>
#include <cstring>
#include <iostream>

//#include "quartet.h"
//#include "locus.h"
//#include "allele.h"
#include "../typedef.h"

/// a class that handels reading/writing the fixed width text format of data.
class format {
public:
	virtual size_t size()=0;			//!< number of characters need for text mode..
	virtual void to_text(uint8_t *, char *)=0;	//!< prints to a row in text mode.
	virtual void from_text(char *, uint8_t *)=0;	//!< reads from a row in text mode.
};


class CHAR7 : public format {
public:
	size_t size(void){return 7;};
	void to_text(uint8_t *, char *);
	void from_text(char *, uint8_t *);
};

class CHAR13 : public format {
public:
	size_t size(void){return 13;};
	void to_text(uint8_t *, char *);
	void from_text(char *, uint8_t *);
};

class UINT8 : public format {
public:
	size_t size(void){return 7;};
	void to_text(uint8_t *, char *);
	void from_text(char *, uint8_t *);
};

#endif
