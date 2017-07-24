#include "sample_name.h"

const std::string Sample_name::file_name=".txt";
const std::string Sample_name::table_name="FILES";
const bool Sample_name::binary=false;

const Registration Sample_name::registered=Registration(Sample_name::table_name, Sample_name::create);

Sample_name::Sample_name ()
{
	mpileup_name="";
	sample_names.clear();
	delim='\t';
}

Sample_name::Sample_name (const std::string &name, const float_t &number)
{
	mpileup_name="";
	sample_names.clear();
	delim='\t';
}

void Sample_name::read (std::istream& in) 
{
	sample_names.clear();
	std::string line;
	std::getline(in, line);
	sample_names=split(line, delim);
	mpileup_name=sample_names.front();
	sample_names.erase(sample_names.begin() );
}

void Sample_name::write(std::ostream& out) const
{
	out << mpileup_name;
	std::vector <std::string>::const_iterator it=sample_names.begin();
	while(it!=sample_names.end() ){
		out << delim << *it;
		++it;
	}
}

std::string Sample_name::header(void) const 
{
	return std::string("@FILNAME\tSMPNAME\t...\n");
}

const std::string Sample_name::sql_header(void) const {
	return "(FILNAME varchar(255), SMPNAME varchar(255) )";
}

const std::string Sample_name::sql_column_names(void) const {
	return "(FILNAME, SMPNAME)";
}

const std::string Sample_name::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
	char *write_ptr=return_buffer;
        std::vector <std::string>::const_iterator it=sample_names.begin();
	//if (it!=sample_names.end() ) write_ptr+=snprintf(return_buffer, SQL_LINE_SIZE, "('%s','%s')", sanitize(mpileup_name).c_str(), sanitize(*it).c_str() );
	if (it!=sample_names.end() ) write_ptr+=snprintf(return_buffer, SQL_LINE_SIZE, "('%s','%s')", mpileup_name.c_str(), it->c_str() );
	++it;
       	while(it!=sample_names.end() ){
		//write_ptr+=snprintf(write_ptr, SQL_LINE_SIZE-(int)(write_ptr-return_buffer), ", ('%s','%s')", sanitize(mpileup_name).c_str(), sanitize(*it).c_str() );
		write_ptr+=snprintf(write_ptr, SQL_LINE_SIZE-(int)(write_ptr-return_buffer), ", ('%s','%s')", mpileup_name.c_str(), it->c_str() );
                ++it;
        }
	return std::string(return_buffer);
}

size_t Sample_name::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+mpileup_name.size()*sizeof(char);
}

const std::string Sample_name::get_file_name(void) const
{
	return file_name;
}

const std::string Sample_name::get_table_name(void) const
{
	return table_name;
}

const bool Sample_name::get_binary(void) const
{
	return binary;
}
