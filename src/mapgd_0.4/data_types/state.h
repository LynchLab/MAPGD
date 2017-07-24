#ifndef _STATE_H_
#define _STATE_H_

#include <string.h>
#include <iostream>
#include <sstream>
#include "bgt.h"
#include "data.h"
#include "typedef.h"

/// A compact containing trait states. Does not store position data.
class State : public data {
private:
	uint32_t *it_, end_;
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new State(Columns);
	}

	size_t char_size_;
	size_t char_per_line_;
	size_t line_length_;
	size_t padding_;
	size_t bit_;
	size_t size_;
	bool bgt_compression_;
	Bgt bgt_;
public:
//Data block
	uint32_t *data;
//end Data block
	State();						//!< simple constructor.
	State(const std::vector <std::string> &);		//!< constructor needed by map_file. String should be column names. 
	State(const State &); 					//!< constructor using a state 
	~State();						//!< destructor.
	size_t size() const;					//!< Returns the number of samples.
	std::string header(void) const;				//!< print header.

	static const std::string table_name;			//!< destination table in Db.
	static const std::string file_name;			//!< default file extension.
	static const bool binary;				//!< flag that indicates whether binary output is avalible.
	get() const;						
	void set(const &,const &);
	void seek(char *);
	char *start(void);
	get_next(void);
	set_next(const uint8_t &);
	set_next8(const uint8_t &);
	const bool get_binary() const;

	State & operator= (const State&);

	void write (std::ostream&) const;
	void read (std::istream&);

	void write_binary (std::ostream&) const;
	void read_binary (std::istream&);
};
#endif
