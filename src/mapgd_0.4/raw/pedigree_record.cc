#include "pedigree_record.h"

//TODO fix these up to make them prettier
std::ostream& operator<< (std::ostream& out, const quartet& q) 
{
	std::stringstream ssout;
	std::string strout;
	out << family << delim << name << delim;

	if (parent[0]!=NULL) out << parent[0]->name << delim;
	else out << delim;

	if (parent[1]!=NULL) out << parent[1]->name << delim;
	else out << delim;
	out << sex;
        return out;
}

std::istream& operator >> (std::istream& in, quartet& q) 
{
	std::string line;
	pedigree->get_record(name);
        return in;
}
