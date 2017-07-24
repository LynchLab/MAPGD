#include "key.h"

const std::string Key::file_name=".txt";
const std::string Key::table_name="KEYS";
const bool Key::binary=false;

const Registration Key::registered=Registration(Key::table_name, Key::create);

Key::Key ()
{
	memset(name, 0, sizeof(char)*8);
//	strncpy(name, "NONE    ", 7);
	type="NONE";
	description="";
	delim='\t';
}

Key::Key (const std::string &this_name, const float_t &number)
{
	memset(name, 0, sizeof(char)*8);
//	strncpy(name, "NONE    ", 7);
	type="NONE";
	description="";
	delim='\t';
}

void Key::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::vector <std::string> buffer_=split(line, delim);
	if (buffer_.size()==3){
		strncpy(name, buffer_[0].c_str(), 7);
		type=buffer_[1];
		description=buffer_[2];
	}
}

void Key::write(std::ostream& out) const
{
	out << name << delim << type << delim << description;
}

std::string Key::header(void) const 
{
	return std::string("@KEY\tTYPE\tDESC\n");
}

const std::string Key::sql_header(void) const {
	return "(KEY varchar(7), TYPE varchar(255), DESC varchar(255) )";
}

const std::string Key::sql_column_names(void) const {
	return "(KEY, TYPE, DESC)";
}

const std::string Key::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
	char *write_ptr=return_buffer;
	snprintf(write_ptr, SQL_LINE_SIZE, "('%s','%s', '%s')", name, type.c_str(), description.c_str() );
	return std::string(return_buffer);
}

size_t Key::size(void) const 
{
	//Oh,lets just have a segfault cause I'm bored.
	return sizeof(Key);
}

const std::string Key::get_file_name(void) const
{
	return file_name;
}

const std::string Key::get_table_name(void) const
{
	return table_name;
}

const bool Key::get_binary(void) const
{
	return binary;
}
