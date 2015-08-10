/* this is a breif demonstration of how I implement ...*/

#include <cstdlib> 
#include <math.h>
#include <initializer_list>
#include <map>
#include <list>
#include <cstring>
#include <iostream>
#include <sqlite3.h>
#include "stream_tools.h"


#define BUFFER 500
#define OMP

class key;
class row;
const void * fetch(const row &,const key &);
void * fetch(row &,const key &);
void * place(row &, key &);
char *error_message = 0;

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


///The key class helps moves data into and out of rows. 
class key
{
private:
	size_t offset_;			//!< specifies where an object is stored relative to the begining of a row.
	size_t size_;			//!< the size of an object in bytes.
	void *begin_;			//!< the reference to where the object is stored.
	bool delete_;			//!< a flag to indicate that key should delete the object at value_ when it goes out of scope.
	bool indexed_;			//!< a flag to indicate that key should be used as an index in the construction of SQL database. 
	std::string name_;		//!< the name of the key.
	std::string desc_;		//!< a verbal description of the key.
	uint8_t keyid_;			//!< the numeric id of the key.
	uint8_t typeid_;		//!< the numeric id of the key.
public:
	/// the default constructor. It basically just zeros everything, and makes sure that we know that the key is empty.

/*	static list <key> keys;
	static list <key> key_names;
	static list <key> key_ids;
*/

	key(){
		name_="";
		begin_=NULL;
		delete_=false;
		size_=0;
		offset_=0;
		typeid_=-1;		//!< the numeric id of the key.
	}

	key(const double &value, const std::string &name, const std::string &desc){
		size_=sizeof(double);
		typeid_=0;
		name_=name;				// set the name of the key.
		desc_=desc;				// set the description of the key.
		begin_=malloc(size_);			// allocate some memory to store the object.
		memcpy(begin_, &value, size_);		// copy over the object to the allocated memory.
		delete_=true;				// toggle the flag that tells us to call free iff the destroctor is called.
		offset_=0;				// the offset needs to be set to some non-sense value.
	};

	/// a constructor that actually makes a meaningful key.
	key(const void *value, const size_t size, const std::string &name, const std::string &desc){
		name_=name;				// set the name of the key.
		desc_=desc;				// set the description of the key.
		begin_=malloc(size);			// allocate some memory to store the object.
		memcpy(begin_, value, size);		// copy over the object to the allocated memory.
		delete_=true;				// toggle the flag that tells us to call free iff the destroctor is called.
		size_=size;				// set the size of the key in bytes.
		offset_=0;				// the offset needs to be set to some non-sense value.
		keyid_=0;
	};

        key(const void *value, const size_t size, const std::string &name){

                name_=name;                             // set the name of the key.
                begin_=malloc(size);                    // allocate some memory to store the object.
                memcpy(begin_, value, size);            // copy over the object to the allocated memory.
                delete_=true;                           // toggle the flag that tells us to call free iff the destroctor is called.
                size_=size;                             // set the size of the key in bytes.
                offset_=0;                              // the offset needs to be set to some non-sense value.
        };

	/// the destructor for key needs to free up the name and description, and the object if it hasn't been freed already.
	~key(){
		clear();
	}

	/// returns the offset from the beging of a row.
	size_t & offset(){
		return offset_;
	}

	/// returns the name of the key.
	std::string & name(){
		return name_;
	}

	/// returns the size of the key in bytes.
	size_t size(void) {
		return size_;
	}

	/// returns a void pointer to the object stored in the key. 
	// this specifically 
	void * value(void) {
		return begin_;
	}

	/// this operator has to be overloaded because we have to keep track of 
	/// whether or not key has memory allocated for an object.
	key & operator=(const key &rhs){
		offset_=rhs.offset_;
		name_=rhs.name_;		
		typeid_=rhs.typeid_;		//!< the numeric id of the key.
		if (rhs.delete_) {		//!< This statement checks to see if rhs has an object 
			if (!delete_) {
				begin_=malloc(rhs.size_);
				delete_=true;
				memcpy(begin_, rhs.begin_, rhs.size_);
			} else {
				if (size_!=rhs.size_) {
					clear();
					begin_=malloc(rhs.size_);
				}	
				memcpy(begin_, rhs.begin_, size_);
			}
		}
		else if (delete_ && !rhs.delete_){
			clear(rhs.begin_);
		}
		size_=rhs.size_;
	}

	/// just a simple check to see if begin_ needs to be freed.
	void clear(void){
		if (delete_){
			free (begin_);
			delete_=false;
			begin_=NULL;
		}
	}

	/// clears begin_, and sets it equal ptr.
	void clear(void *ptr){
		clear();
		begin_=ptr;
	}

	/// the copy construcotr. The only thing that needs to be done is to tell 
	/// the key that it doesn't have any memory allocated to it yet.
	key(key const& rhs) {
		delete_=false; 
		*this=rhs;
	}

	///The basic way we move info ... 
	const char * to_string(const void *ptr) const{
		switch (typeid_){
			case (0):
				return std::to_string(*(double *)(ptr) ).c_str();
			case (1):
				return std::to_string(*(int *)(ptr) ).c_str();
			//case (3):
			//	return std::to_string( (char *)(ptr) );
			default:
				break;
		}
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed to print object.\n";
		exit(0);
		return "";
	};
	

	friend const void * fetch(const row &, const key &);
	friend void * fetch(row &, const key &);
	friend void * place(row &, key &);
};

///The row class...
class row
{
private:
	char *begin_;
	size_t size_;
	std::string table_name_;
	std::list <key> keys_;
	std::map <const std::string, key *> key_map;
public:
	row(std::string table_name, std::initializer_list <key> this_list){
		table_name_=table_name;
		size_=0;
		begin_=(char *)malloc(size_);
		for(key this_key : this_list ) place(*this, this_key);
	};
	~row(void){
		free(begin_);
	}
	key get_key(const std::string &name) const {
		std::map<std::string, key*>::const_iterator it = key_map.find(name);
		return *(it->second );
	}

	friend const void * fetch(const row &, const key &);
	friend void * fetch(row &,const key &);
	friend void * place(row &, key &);
	friend void * db_insert(sqlite3 *, const row &);

	const char * table(void) const{
		return table_name_.c_str();
	}
};


inline void * fetch(row &this_row, const key &this_key)
{
	return this_row.begin_+this_key.offset_;
}

inline const void * fetch(const row &this_row, const key &this_key) 
{
	return this_row.begin_+this_key.offset_;
}

void * place(row &this_row, key &this_key)
{
	key *row_key=this_row.key_map[this_key.name()];
	if( row_key == NULL )
	{
		this_row.keys_.push_back(this_key);
		this_row.key_map[this_key.name()]=&(this_row.keys_.back());
		row_key=&(this_row.keys_.back() );
	}
	char *new_block=(char *)malloc(this_row.size_+row_key->size_);
	if (new_block==NULL) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed place object in row, try decreasing your buffer size and recompiling.\n";
	}
	memcpy(new_block, this_row.begin_, this_row.size_);
	memcpy(new_block+this_row.size_, row_key->begin_, row_key->size_);
	row_key->offset_=this_row.size_;
	this_key.clear(this_row.begin_+this_key.offset_);
	free(this_row.begin_);
	this_row.begin_=new_block;
	this_row.size_=this_row.size_+row_key->size_;
	return this_row.begin_+this_key.offset_;
}

/*
void * db_fetch(const id0_t &chrom, const id1_t &bp, const key &chrom, const key &get_me, sqlite3 &db)
{
	std::string sql=strprintf("SELECT %s FROM GENOME WHERE ROWID=%d AND CHROM=%d); ", get_me.name(), bp, chrom);
	sqlite3_exec(db, sql.c_str(), callback, 0, &error_mesage);
}


void * db_update(const id0_t &chrom, const id1_t &bp, sqlite3 &db)
{
	std::string sql=strprintf("SELECT %s FROM GENOME WHERE ROWID=%d AND CHROM=%d); ", get_me.name(), bp, chrom);
	sqlite3_exec(db, sql.c_str(), callback, 0, &error_mesage);
}*/

void * db_insert(sqlite3 *db, const row &this_row)
{
	std::string make_table=strprintf("CREATE TABLE IF NOT EXIST %s ", this_row.table() );
	std::string sql=strprintf("INSERT INTO %s VALUES (", this_row.table() );
	std::cout << sql << std::endl;
	for ( auto &this_key : this_row.keys_){
		if (&this_key!=&this_row.keys_.back() ) sql=strprintf("%s %s,", sql.c_str(), this_key.to_string(fetch(this_row, this_key) ) );
		else sql=strprintf("%s %s", sql.c_str(), this_key.to_string(fetch(this_row, this_key) ) );
		std::cout << sql << std::endl;
	};
	sql=strprintf("%s);", sql.c_str() );
	std::cout << sql << std::endl;
	sqlite3_exec(db, sql.c_str(), callback, 0, &error_message);
}


int main (int argc, char *argv[])
{
	sqlite3 *db;
	int rc=sqlite3_open("test.db", &db);

	double x=0, y=1;
	key result(x, std::string("RES"), std::string("The sign of the log") );

	key value(y, std::string("ONE"), std::string("The number one") ); 

	double *res, *val; 

	const row input[BUFFER]=row("input", {value});
	row output[BUFFER]=row("output", {result});

//	db_insert(db, input[0]);
//	db_insert(db, output[0]);

	result=output[0].get_key("RES");
	value=input[0].get_key("ONE");

std::string   sql = "CREATE TABLE GENOME("  \
         "ROW INT PRIMARY KEY     NOT NULL," \
         "MAJOR          TEXT     NOT NULL," \
         "MINOR          TEXT     NOT NULL," \
         "FREQ           REAL     NOT NULL," \
         "ERROR          REAL     NOT NULL," \
         "LOGLIKE        REAL     NOT NULL);";

   /* Execute SQL statement */
   rc = sqlite3_exec(db, sql.c_str(), callback, 0, &error_message);


	for (int y=0; y<50000; ++y){
#ifdef OMP
		#pragma omp parallel private(res, val) 
		#pragma omp for
#endif
		for (int x=0; x<BUFFER; ++x){
			res=(double *)fetch(input[x], result);
			val=(double *)fetch(output[x], value);
			*res=sin( log (*val) );
		}
	}
	sqlite3_close(db);
	std::cout << "done.\n";
}
