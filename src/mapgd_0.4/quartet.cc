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

std::ostream& operator<< (std::ostream& out, const quartet& q) {
	if (not(q.masked) ) out << int(q.base[0]) << q.delim << int(q.base[1]) << q.delim << int(q.base[2]) << q.delim << int(q.base[3]); 
	else out << "0" << q.delim << "0" << q.delim << "0" << q.delim << "0"; 
        return out;
};

