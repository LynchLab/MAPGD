#include "quartet.h"

count_t count(const quartet_t &q){
	return q.base[0]+q.base[1]+q.base[2]+q.base[3];
}

void unmask(quartet_t &q){
	q.masked=false;
}

void mask(quartet_t &q){
	q.masked=true;
}

//TODO fix these up to make them prettier
std::ostream& operator<< (std::ostream& out, const quartet& q) 
{
	std::stringstream ssout;
	std::string strout;
	if (not(q.masked) ){ 
		ssout << int(q.base[0]);
		ssout << q.delim;
		ssout << int(q.base[1]);
		ssout << q.delim;
		ssout << int(q.base[2]);
		ssout << q.delim;
		ssout << int(q.base[3]);
		strout=ssout.str();
		out << strout;
		if (strout.size()<15)
			out << std::string(15-strout.size(), ' ');
	}	
	else out << "0" << q.delim << "0" << q.delim << "0" << q.delim << "0     "; 
        return out;
}

//TODO fix these up to make them prettier
std::istream& operator >> (std::istream& in, quartet& q) 
{
	std::string line;
	getline(in, line, '/');
	q.base[0]=atoi(line.c_str() );
	getline(in, line, '/');
	q.base[1]=atoi(line.c_str() );
	getline(in, line, '/');
	q.base[2]=atoi(line.c_str() );
	getline(in, line, '\t');
	q.base[3]=atoi(line.c_str() );
        return in;
}
