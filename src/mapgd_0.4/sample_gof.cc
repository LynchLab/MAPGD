#include "sample_gof.h"

const std::string Sample_gof::file_name=".gof";
const std::string Sample_gof::table_name="SAMPLE";

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

std::istream& operator >> (std::istream& in, Sample_gof& x)
{
	in >> x.name_ >> x.number_;
	return in;
}

std::ostream& operator<< (std::ostream& out, const Sample_gof& x) 
{
	out << x.name_ << x.delim << x.number_;
	return out;
}

std::string Sample_gof::header(void) const 
{
	return "@SMPNAME\tGOF\n";
}

size_t Sample_gof::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+name_.size()*sizeof(char);
}
