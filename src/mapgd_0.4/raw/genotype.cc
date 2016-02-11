#include "genotype.h"
Genotype::Genotype()
{
	lMM=-REAL_MAX; 
	lMm=-REAL_MAX; 
	lmm=-REAL_MAX; 
	N=0;
}
Genotype::Genotype(const float_t &MM, const float_t &Mm, const float_t &mm, const count_t &lN)
{
	lMM=MM; 
	lMm=Mm; 
	lmm=mm; 
	N=lN;
}
	
Genotype & Genotype::operator= (const Genotype& rhs)
{
	lMM=rhs.lMM;
	lMm=rhs.lMm;
	lmm=rhs.lmm;						//!< Major Major, Major minor, minor minor
	N=rhs.N;						//!< total depth of coverage.
	return *this;
}

std::ostream& operator<< (std::ostream& out, const Genotype& x)
{
	out << x.lMM << '/' << x.lMm << '/' << x.lmm << '/' << x.N;
	return out;
}

std::istream& operator>> (std::istream& in, Genotype& x)
{
	std::string line;
	getline(in, line, '/');
	x.lMM=std::stof(line);
	getline(in, line, '/');
	x.lMm=std::stof(line);
	getline(in, line, '/');
	x.lmm=std::stof(line);
	getline(in, line, '\t');
	x.N=std::stoi(line);
        return in;
}
