#ifndef _EXTERNAL_DATA_H_
#define _EXTERNAL_DATA_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include <stdarg.h>

#include "typedef.h"

//! An interface used to read/write data from outside of mapgd.
/*! This interface has to associate external data with internal data
    so that mapgd only has to touch the internal data classes. 
 */
class External_data {

public:
	//! The write function must be defined in the child class.
	virtual void put(const Data* data, ...) = 0;
	//! The read function must be defined in the child class.
	virtual void get(Data* data, ...) const = 0;

};

#endif
