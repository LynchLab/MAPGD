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

std::ostream& operator<<(std::ostream& output, const quartet& q)
{
	output << int(q.base[0]) << q.delim << int(q.base[1]) << q.delim << int(q.base[2]) << q.delim << int(q.base[3]);
	return output;
}

std::istream& operator>>(std::istream& input, quartet& q)
{
	int temp;
	input >> temp;
	q.base[0]=count_t(temp);
	input.ignore(1, q.delim);
	
	input >> temp;
	q.base[1]=count_t(temp);
	input.ignore(1, q.delim);

	input >> temp;
	q.base[2]=count_t(temp);
	input.ignore(1, q.delim);

	input >> temp;
	q.base[3]=count_t(temp);
	return input;
}
