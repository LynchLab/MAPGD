#ifndef _FORMAT_THIS_H_
#define _FORMAT_THIS_H_

#include <string>
#include <cstring>
#include <iostream>

//#include "quartet.h"
//#include "locus.h"
//#include "allele.h"
#include "../typedef.h"

/// a class that handels reading/writing the fixed width text format of data.
// Since we don't know what kind of data is referenced from the base class, 
// there is nothing we can do and the function just advances the stream

void to_text_char7 (uint8_t * txt, void* val); 
void from_tx_char7 (void *val, uint8_t* txt);
size_t size_of_char7 (void *);

void to_text_char13 (uint8_t * txt, void* val); 
void from_tx_char13 (void *val, uint8_t* txt);
size_t size_of_char13 (void *);

void to_text_uint8 (uint8_t * txt, void* val);
void from_tx_uint8 (void *val, uint8_t* txt);
size_t size_of_uint8 (void *);

size_t size_of_real (void *);
void to_text_real (uint8_t * txt, void* val);
void from_tx_real (void *val, uint8_t* txt);
	
size_t size_of_null (void *);
void to_text_null (uint8_t*, void*);
void from_tx_null (void*, uint8_t*);

#endif
