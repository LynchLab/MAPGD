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

protected:
	//! The read function must be defined in the child class.
	virtual void read(std::istream&  , Data& data, ...) = 0;
	//! The write function must be defined in the child class.
	virtual void write(std::ostream& , Data& data, ...) const = 0;
public:
	External_data(){};
	External_data(std::vector <std::string> &){};

	/**
	 * \defgroup basic IO
	 * @{
	 */
	//! Use the \<\< operator to write Data in text mode.
	friend std::ostream& operator << (std::ostream&, const Data &);	
	//! Use the \>\> operator to read Data in text mode.
	friend std::istream& operator >> (std::istream&, Data &);		

	//! The size of the class in bytes. 
	/* If size() is undefined in the child class, then we don't want to 
	 * read or write any at all, so we return 0.
	 */
        virtual size_t size(void) const {return 0;};
	/*
	 *@}
         */

};

#endif
