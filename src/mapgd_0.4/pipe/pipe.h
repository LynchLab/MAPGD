/* this is a breif demonstration of how I implement ...*/

#ifndef _PIPE_H_
#define _PIPE_H_

#include <cstdlib> 
#include <math.h>
#include <initializer_list>
#include <map>
#include <list>
#include <cstring>
#include <iostream>
#include <sqlite3.h>
#include "../typedef.h"

class key;
class row;
class table;

void * fetch(const row &, const key &);
void * place(row &, const key &, void *);
void * reserve(row &, const key &);

int write_row(std::ostream &, const row &);
int read_row(std::istream&, row &);


/// the key class lets us move data into and out of rows. 
/* 	Each key has a breif description of the data it represents, a ...
 */
class key
{
private:
	size_t offset_;			//!< specifies where an object is stored relative to the begining of a row.
	size_t size_;			//!< the size of an object in bytes.
	bool indexed_;			//!< a flag to indicate that key should be used as an index in the construction of a database. 
	char name_[256];		//!< the name of the key.
	char desc_[256];		//!< a verbal description of the key.
	keyid_t keyid_;			//!< the numeric id of the key.
	keyid_t typeid_;		//!< the numeric id of the data type. May make this type_info.
public:
	static const keyid_t nokey=-1;
	/// the default constructor. It basically just zeros everything, and makes sure that we know that the key is empty.
	key();
	
	/// overloading
	key(const keyid_t &type_id, const char *name, const char *desc);

	/// a constructor that actually makes a meaningful key.
	key(const keyid_t &type_id, const size_t &size, const char *name, const char *desc);
	/// returns the offset from the beging of a row.
	size_t & offset();
	/// returns the name of the key.
	const char * name() const;
	/// returns the size of the key in bytes.
	size_t size(void);
	/// returns a void pointer to the object stored in the key. 
	// this specifically 
	void * value(void);
	/// this operator has to be overloaded because we have to keep track of 
	/// whether or not key has memory allocated for an object.
	key & operator=(const key &rhs);
	/// returns a textual representation of the data
	const char * to_string(const void *ptr) const;
	/// returns the SQL type of a key.
	const char * sql_type (void) const;
	
	friend void * fetch(const row &, const key &);
	friend void * place(row &, const key &, void *);
	friend void * reserve(row &, const key &);
};

/// The row class stores the data for inout/output operations.
/* 	Right now each row stores a list of keys, this will probably be moved over to a table class.
 *
 */
class row
{
private:
	char *begin_;					//!< where the memory is
	size_t size_;					//!< amount of memory allocated
	std::string table_name_;			//!< name of the table to which the row belongs
	std::list <key> keys_;				//!< list of keys (none of which have size set to 0)
	std::map <const std::string, key *> key_map;	//!< a hash table of key names.
public:
	row(std::string table_name, std::initializer_list <key> this_list);	//!< initilizes a row with a list of keys.
	~row(void);								//!< frees up some memory.
	key get_key(const std::string &name) const;				//!< returns a key from the hash table.
	friend void * fetch(const row &, const key &);
	friend void * place(row &, const key &, void *);
	friend void * reserve(row &, const key &);

	friend int read_row(std::istream &, row &);
	friend int write_row(std::ostream &, const row &);

	friend void * db_insert(sqlite3 *, const row &);
	const char * table(void) const;
};

class table
{
private:
	std::string table_name_;			//!< name of the table to which the row belongs
	std::list <key> keys_;				//!< list of keys (none of which have size set to 0)
	std::map <const std::string, key *> key_map;	//!< a hash table of key names.
public:
	table(std::string table_name, sqlite3 *);			//!< initilizes a row with a list of keys.
	key get_key(const std::string &name) const;			//!< returns a key from the hash table.
	void * write(sqlite3 *, const row &);
};

inline void * fetch(const row &this_row, const key &this_key)
{
	return this_row.begin_+this_key.offset_;
}

inline int write_row(std::ostream &out, const row &this_row)
{
	if (out.write(this_row.begin_, this_row.size_ ) ) return 1; 
	else return 0;
}

inline int read_row(std::istream &in, row &this_row)
{
	if (in.read(this_row.begin_, this_row.size_ ) ) return 1; 
	else return 0;
}

#endif
