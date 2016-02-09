#include "genotype.h"
#include <iostream>
using namespace std;

const std::string population_genotypes::file_name="gcf";
const std::string population_genotypes::table_name="GENOTYPES";

/** @breif default constructor. Does nothing **/

genotype::genotype(){}
/** @breif constuctor w/ initial values. **/

genotype::genotype(const float_t &MM, const float_t &Mm, const float_t &mm, const count_t &lN)
{
	lMM=MM; 
	lMm=Mm; 
	lmm=mm; 
	N=lN;
}
	
population_genotypes::population_genotypes(const std::vector <std::string> &column_names)
{
	sample_names_=std::vector <std::string> (column_names.cbegin()+3, column_names.cend() );
	likelihoods.resize(sample_names_.size() );
}

population_genotypes::population_genotypes()
{
	frozen=false;
}

population_genotypes::population_genotypes(const population_genotypes &rhs)
{
	*this=rhs;
}

population_genotypes::~population_genotypes(){}

/**@breif return size of population_genotypes if population_genotypes is set, 0 otherwise**/
size_t population_genotypes::size() const
{
	if (frozen) return likelihoods.size();
	else return 0;
}

void population_genotypes::add(const genotype &_gl)
{
	if (frozen){
		*igl_=_gl;
		igl_++;
	} else {
		likelihoods.push_back(_gl);	
	};
}

void population_genotypes::add(const float_t &lMM, const float_t &lMm, const float_t &lmm, const count_t &lN)
{
	if (frozen){
		*igl_=genotype(lMM, lMm, lmm, lN);
		igl_++;
	} else {
		likelihoods.push_back(genotype(lMM, lMm, lmm, lN));	
	};
}

void population_genotypes::clear()
{
	frozen=true;
	igl_=likelihoods.begin();
}

Genotype_pair_tuple convert(const genotype &x, const genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
	return Genotype_pair_tuple ( roundf(x.lMM * T) / T, roundf(x.lMm*T)/T, roundf(x.lmm*T)/T, roundf(y.lMM*T)/T, roundf(y.lMm*T)/T, round(y.lmm*T)/T, round(m*T)/T);
}

std::string population_genotypes::header(void) const
{
	std::string line="@SCFNAME\tPOS\tMN_FREQ";
	std::vector <std::string>::const_iterator s_it=sample_names_.cbegin(), end=sample_names_.cend();
	while(s_it!=end){
		line+='\t';	
		line+=*s_it;
		s_it++;	
	}
	line+='\n';
	return line;
}

genotype & genotype::operator= (const genotype& rhs)
{
	lMM=rhs.lMM;
	lMm=rhs.lMm;
	lmm=rhs.lmm;						//!< Major Major, Major minor, minor minor
	N=rhs.N;						//!< total depth of coverage.
	return *this;
}

std::ostream& operator<< (std::ostream& out, const population_genotypes& x)
{
	out << x.id1 << '\t' << x.m;
	std::vector <genotype>::const_iterator s_it=x.likelihoods.cbegin(), end=x.likelihoods.cend();
	while(s_it!=end){
		out << '\t' << *s_it;
		s_it++;	
	}
	return out;
}

std::istream& operator >> (std::istream& in, population_genotypes& x)
{
        std::string line;
        std::getline(in, line);
        std::stringstream line_stream(line);

	line_stream >> x.id1 >> x.m;
	std::vector <genotype>::iterator s_it=x.likelihoods.begin(), end=x.likelihoods.end();
	while(s_it!=end){
		line_stream >> *s_it;
		s_it++;
	}	
	return in;
}

std::ostream& operator<< (std::ostream& out, const genotype& x)
{
	out << x.lMM << '/' << x.lMm << '/' << x.lmm << '/' << x.N;
	return out;
}

std::istream& operator>> (std::istream& in, genotype& x)
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

Genotype_pair::Genotype_pair(const float_t &X_MM_, const float_t &X_Mm_, const float_t &X_mm_, const float_t &Y_MM_, const float_t &Y_Mm_, const float_t &Y_mm_, const float_t &m_)
{
	X_MM=X_MM_;
	X_Mm=X_Mm_;
	X_mm=X_mm_;
	Y_MM=Y_MM_;
	Y_Mm=Y_Mm_;
	Y_mm=Y_mm_;
	m=m_;

}

Genotype_pair_tuple Genotype_pair::to_tuple(const Genotype_pair &pair)
{
	return Genotype_pair_tuple (pair.X_MM, pair.X_Mm, pair.X_mm, pair.Y_MM, pair.Y_Mm, pair.Y_mm, pair.m);
}

Genotype_pair Genotype_pair::from_tuple(const Genotype_pair_tuple &t)
{
	return Genotype_pair (std::get<0>(t),std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), std::get<5>(t), std::get<6>(t) );
}

population_genotypes & population_genotypes::operator= (const population_genotypes& rhs)
{
	sample_names_=rhs.sample_names_;	//!< a vector of sample names.
	frozen=rhs.frozen;					//!< a flag to indicate that no more samples will be added 

	likelihoods=rhs.likelihoods;		//!< genotypic likelihood
	igl_=likelihoods.begin();			//!< an iterator to allow us to iterate over the likelihoods.

	major=rhs.major;				//!< identity of the major allele
	minor=rhs.minor;				//!< identity of the minor allele
	m=rhs.m;					//!< minor allele frequency
	id0=rhs.id0;					//!< scaffold number
	id1=rhs.id0;					//!< location on scaffold
	return *this;
}
