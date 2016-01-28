#include "clone_gof.h"

const std::string Clone_gof::file_name=".gof";
const std::string Clone_gof::table_name="SAMPLE";

Clone_gof::Clone_gof (){
	name_="";
	number_=0;
	delim='\t';
}

Clone_gof::Clone_gof (const std::string &name, const float_t &number){
	name_=name;
	number_=number;
	delim='\t';
}

std::istream& operator >> (std::istream& in, Clone_gof& x) {
	in >> x.name_ >> x.number_;
	return in;
}

std::ostream& operator<< (std::ostream& out, const Clone_gof& x) {
	out << x.name_ << x.delim << x.number_;
	return out;
};

std::string Clone_gof::header(void) const {
	return "@SAMPLE\tGOF\n";
}

size_t Clone_gof::size(void) const {
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	
	float_t number_;
	std::string name_;
	
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+name_.size()*sizeof(char);
}
