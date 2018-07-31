/* An bed  */

#ifndef _BED_FILE_H_
#define _BED_FILE_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#ifndef NOLZ4
#include "state.h"
#endif

#include "external-file.h"
#include "external-data.h"

#include "typedef.h"
#include "datatypes.h"
#include "raw.h"
#include "stream_tools.h"

//Bed files do not contain sufficient information to be interpreted without other files (because they are terribly designed,

class Bed_data : public External_data {
private:
public:
	size_t size_;
	uint8_t *buffer_;
	Bed_data ();
	Bed_data (const size_t &);

	void put (const Data *, ...);
	void get (Data *, ...) const;

#ifndef NOLZ4
	id1_t get (State &) const;
#endif 

};


class Bed_file : public External_file <Bed_data> {
private:
	using Base_file::in_;
	using Base_file::out_;
	using Base_file::open_;
	using Base_file::open_no_extention;
//	htsFile *file_;
public:
	void open(const std::ios_base::openmode &);
	void open(const char *, const std::ios_base::openmode &);
	void close(void);
	void read (Bed_data &);
};

#endif
