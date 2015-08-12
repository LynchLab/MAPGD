/* this is a breif demonstration of how I implement ...*/

#include "pipe.h"

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


/// the default constructor. It basically just zeros everything, and makes sure that we know that the key is empty.
key::key(){
	memset(name_, 0, 255);	// this key has no name.
	memset(desc_, 0, 255);	// this key has no name.
	begin_=NULL;		// the data this key refernces is stored at NULL.
	delete_=false;		// if this key is deconstructed, then it won't free any memory.
	size_=0;		// if this key is called, then 0 bytes will be read from the data stream.
	offset_=0;		// the data for this key stats 0 bytes after the begging of a row.
	keyid_=nokey;
	typeid_=nokey;		// a nonsense type id.
}

//TODO write some python code to update this function.
key::key(const keyid_t &typeid_, const char *name, const char *desc){
	switch (typeid_){
		case key_type::real_32:	key(key_type::real_32, 32, name, desc);
		case key_type::real_64:	key(key_type::real_64, 64, name, desc);
		case key_type::uint_8:	key(key_type::uint_8,   8, name, desc);
		case key_type::uint_16:	key(key_type::uint_16, 16, name, desc);
		case key_type::uint_32:	key(key_type::uint_32, 32, name, desc);
		case key_type::locus: 	key(key_type::locus,   0, name, desc);
		case key_type::allele: 	key(key_type::allele,  0, name, desc);
		case key_type::genotyp:	key(key_type::genotyp, 0, name, desc);
	}
}

/// a constructor that actually makes a meaningful key.
key::key(const keyid_t &type_id, const size_t &size, const char *name, const char *desc){

	strncpy(name_, name, 255);
	strncpy(desc_, desc, 255);
	name_[255]=0;
	desc_[255]=0;

	size_=size;				// set the size of the key in bytes.
	offset_=0;				// the offset needs to be set to some non-sense value.
	keyid_=0;				// we don't know the key id.
	typeid_=type_id;			// they told us the type id.
}

/// returns the offset from the beging of a row.
size_t & key::offset(){
	return offset_;
}

/// returns the name of the key.
const char * key::name() const{
	return name_;
}

/// returns the size of the key in bytes.
size_t key::size(void) {
	return size_;
}

/// this operator has to be overloaded because we have to keep track of 
/// whether or not key has memory allocated for an object.
key & key::operator=(const key &rhs){
	offset_=rhs.offset_;
	memcpy(name_, rhs.name_, 256);	
	typeid_=rhs.typeid_;		//!< the numeric id of the key.
	size_=rhs.size_;
	return *this;
}

/// returns a textual representation of the data
const char * key::to_string(const void *ptr) const{
	switch (0){
		case key_type::real_32:	std::to_string( *(float *) ptr );
		case key_type::real_64:	std::to_string( *(double *) ptr );
		case key_type::real_80:	std::to_string( *(long double *) ptr );
		case key_type::uint_8:	std::to_string( *(uint8_t *) ptr );
		case key_type::uint_16:	std::to_string( *(uint16_t *) ptr );
		case key_type::uint_32:	std::to_string( *(uint32_t *) ptr );
		case key_type::uint_64:	std::to_string( *(uint64_t *) ptr );
		default:
			std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed to format object.\n";
			break;
	}
	std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed to format object.\n";
	exit(0);
	return "";
};

/// returns the SQL type of a key.
const char * key::sql_type (void) const {
	switch(typeid_){
		case key_type::real_32:	return "REAL (32)";
		case key_type::real_64:	return "REAL (64)";
		case key_type::uint_8:	return "INT (8)";
		case key_type::uint_16:	return "INT (16)";
		case key_type::uint_32:	return "INT (32)";
		case key_type::uint_64:	return "INT (64)";
		case key_type::locus: 	return "";
		case key_type::allele: 	return "";
		case key_type::genotyp:	return "";
	}
	return "";
}


///The row class...
row::row(std::string table_name, std::initializer_list <key> this_list){
	table_name_=table_name;
	size_=0;
	begin_=(char *)malloc(size_);
	for(key this_key : this_list ) reserve(*this, this_key);
}

row::~row(void){
	free(begin_);
}

key row::get_key(const std::string &name) const {
	std::map<std::string, key*>::const_iterator it = key_map.find(name);
	return *(it->second );
}

const char * row::table(void) const{
	return table_name_.c_str();
};


void * place(row &this_row, const key &this_key, void * this_value)
{
	key *row_key=this_row.key_map[this_key.name()];
	if( row_key == NULL )
	{
		this_row.keys_.push_back(this_key);
		this_row.key_map[this_key.name()]=&(this_row.keys_.back());
		row_key=&(this_row.keys_.back() );
	}
	char *new_block=(char *) malloc (this_row.size_+row_key->size_);
	if (new_block==NULL) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed place object in row, try decreasing your buffer size and recompiling.\n";
	}
	memcpy(new_block, this_row.begin_, this_row.size_);
	//TODO insert a check here to see if this_value is of the right type and size.	
	memcpy(new_block+this_row.size_, this_value, row_key->size_);
	row_key->offset_=this_row.size_;
	this_row.begin_=new_block;
	this_row.size_=this_row.size_+row_key->size_;

	return this_row.begin_+this_key.offset_;
}

void * reserve(row &this_row, const key &this_key)
{
	key *row_key=this_row.key_map[this_key.name()];
	if( row_key == NULL )
	{
		this_row.keys_.push_back(this_key);
		this_row.key_map[this_key.name()]=&(this_row.keys_.back());
		row_key=&(this_row.keys_.back() );
	}
	char *new_block=(char *) malloc (this_row.size_+row_key->size_);
	if (new_block==NULL) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: failed place object in row, try decreasing your buffer size and recompiling.\n";
	}
	memcpy(new_block, this_row.begin_, this_row.size_);
	
	memset(new_block+this_row.size_, 0, row_key->size_);
	row_key->offset_=this_row.size_;
	this_row.begin_=new_block;
	this_row.size_=this_row.size_+row_key->size_;
	return this_row.begin_+this_key.offset_;
}
