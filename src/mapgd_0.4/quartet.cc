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
	if (not(q.masked) ){ 
		out << std::setfill('0') << std::setw(3) << int(q.base[0]);
		out << q.delim;
		out << std::setfill('0') << std::setw(3) << int(q.base[1]);
		out << q.delim;
		out << std::setfill('0') << std::setw(3) << int(q.base[2]);
		out << q.delim;
		out << std::setfill('0') << std::setw(3) << int(q.base[3]); 
	}
	else out << "000" << q.delim << "000" << q.delim << "000" << q.delim << "000"; 
        return out;
};

std::istream& operator >> (std::istream& in, quartet& q) {
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
};
