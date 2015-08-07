#include "row.h"

row::row (std::list <key> keys )
{
	std::list <key>::iterator it=keys.begin();
	std::list <key>::iterator end=keys.end();
	row_size_=0;
	while (it!=end){
		key this_key=*it;
		this_key.set_offset(row_size_);
		row_size_+=it->size();
		keys_.push_back(this_key);
		keys_by_name_[it->get_name()]=&keys_.back();
		keys_by_num_[it->get_num()]=&keys_.back();
		++it;
	}
	data_=(uint8_t *)calloc(1, row_size_);
}

row::row ()
{
        row_size_=0;
        data_=(uint8_t *)calloc(1, row_size_);
}

row::~row(void)
{
	free(data_);
}	

void row::add_key (const key &this_key)
{
	if (this_key.is!=key::nokey);	// checking to see if the key has been added to a row.
	size_t old_row_size=row_size_;
	key new_key=this_key;
	new_key.set_offset(row_size_);
	row_size_+=new_key.size();
	keys_.push_back(new_key);
	keys_by_name_[new_key.get_name()]=&keys_.back();
	keys_by_num_[new_key.get_num()]=&keys_.back();
	uint8_t * new_data=(uint8_t *)calloc(1, row_size_);
	if(new_data==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: out of memory.\n";
		exit(0);
	}
	memcpy(new_data, data_, old_row_size);
	memcpy(new_data, &this_key.value(), new_key.size() );
	free(data_);
	this_key.clear();
	data_=new_data;	
}

std::list <key> row::get_keys(void) const
{
	return keys_;
}		
	
key row::get_key(const char *key_name) const
{
	return get_key(std::string(key_name) );
}

key row::get_key(const std::string &key_name) const
{
	return *keys_by_name_.find(key_name)->second; 
}			

key row::get_key(const uint8_t &key_num) const
{
	return *keys_by_num_[key_num]; 
}			

void row::get(const std::string &key_name, void *dest) const
{
	key this_key=get_key(key_name);
	memcpy(dest, data_+this_key.offset(), this_key.size() );
}

void row::get(const uint8_t &key_num, void *dest) const
{
	key this_key=get_key(key_num);
	memcpy(dest, data_+this_key.offset(), this_key.size() );
}

void row::get(const char *key_name, void *dest) const
{
	key this_key=get_key(key_name);
	memcpy(dest, data_+this_key.offset(), this_key.size() );
}
