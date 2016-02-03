#include "sample_name.h"

const std::string Sample_name::file_name=".txt";
const std::string Sample_name::table_name="FILES";

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

std::istream& operator >> (std::istream& in, Sample_name& x)
{
	x.sample_names.clear();
	std::string line;
	std::getline(in, line);
	x.sample_names=split(line, x.delim);
	x.mpileup_name=x.sample_names.front();
	x.sample_names.erase(x.sample_names.begin() );
	return in;
}

std::ostream& operator<< (std::ostream& out, const Sample_name& x) 
{
	out << x.mpileup_name;
	std::vector <std::string>::const_iterator it=x.sample_names.begin();
	while(it!=x.sample_names.end() ){
		out << x.delim << *it;
		++it;
	}
	return out;
}

std::string Sample_name::header(void) const 
{
	return "@FILNAME\tSMPNAME\t...\n";
}

size_t Sample_name::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+mpileup_name.size()*sizeof(char);
}
