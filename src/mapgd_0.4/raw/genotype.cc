#include "genotype.h"
Genotype::Genotype()
{
	MM=0; 
	Mm=0; 
	mm=0; 
	N=0;
}
Genotype::Genotype(const float_t &MM_, const float_t &Mm_, const float_t &mm_, const count_t &N_)
{
	MM=MM_; 
	Mm=Mm_; 
	mm=mm_; 
	N=N_;
}
	
Genotype & Genotype::operator= (const Genotype& rhs)
{
	MM=rhs.MM;
	Mm=rhs.Mm;
	mm=rhs.mm;						//!< Major Major, Major minor, minor minor
	N=rhs.N;						//!< total depth of coverage.
	return *this;
}

std::ostream& operator<< (std::ostream& out, const Genotype& x)
{
	out << x.MM << '/' << x.Mm << '/' << x.mm << '/' << x.N;
	return out;
}

std::istream& operator>> (std::istream& in, Genotype& x)
{
	std::string line;
	getline(in, line, '/');
	//TODO make this respect the type of x.MM
	x.MM=std::stod(line);
	getline(in, line, '/');
	x.Mm=std::stod(line);
	getline(in, line, '/');
	x.mm=std::stod(line);
	getline(in, line, '\t');
	x.N=std::stoi(line);
        return in;
}
