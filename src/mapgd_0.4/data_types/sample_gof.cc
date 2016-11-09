#include "sample_gof.h"

const std::string Sample_gof::file_name=".gof";
const std::string Sample_gof::table_name="SAMPLE";
const bool Sample_gof::binary=false;

const Registration Sample_gof::registered=Registration(Sample_gof::table_name, Sample_gof::create);

Sample_gof::Sample_gof ()
{
	name_="";
	number_=0;
	delim='\t';
}

Sample_gof::Sample_gof (const std::string &name, const float_t &number)
{
	name_=name;
	number_=number;
	delim='\t';
}

void
Sample_gof::read (std::istream& in)
{
	in >> name_ >> number_;
}

void
Sample_gof::write (std::ostream& out) const
{
	out << name_ << delim << number_;
}

std::string 
Sample_gof::header(void) const 
{
	return "@SMPNAME\tGOF\n";
}

size_t 
Sample_gof::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+name_.size()*sizeof(char);
}

const bool Sample_gof::get_binary(void) const
{
	return binary;
}
