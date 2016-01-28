#include "genotype.h"
#include <iostream>
using namespace std;

const std::string population_genotypes::file_name="gcf";
const std::string population_genotypes::table_name="GENOTYPES";

/** @breif default constructor. Does nothing **/

genotype::genotype(){};
/** @breif constuctor w/ initial values. **/

genotype::genotype(const float_t &MM, const float_t &Mm, const float_t &mm, const count_t &lN){
	lMM=MM; 
	lMm=Mm; 
	lmm=mm; 
	N=lN;
}
	
population_genotypes::population_genotypes(std::vector <std::string> column_names)
{
	sample_names=column_names;
	likelihoods.resize(sample_names.size() );
}

population_genotypes::population_genotypes(){frozen=false;};
population_genotypes::population_genotypes(const population_genotypes &poplikelihoods){
	likelihoods=poplikelihoods.likelihoods;
	m=poplikelihoods.m;
	igl=poplikelihoods.igl;
	frozen=poplikelihoods.frozen;
}

population_genotypes::~population_genotypes(){};

/**@breif return size of population_genotypes if population_genotypes is set, 0 otherwise**/
size_t population_genotypes::size() const{
	if (frozen) return likelihoods.size();
	else return 0;
}

void population_genotypes::add(const genotype &_gl){
	if (frozen){
		*igl=_gl;
		igl++;
	} else {
		likelihoods.push_back(_gl);	
	};
}

void population_genotypes::add(const float_t &lMM, const float_t &lMm, const float_t &lmm, const count_t &lN){
	if (frozen){
		*igl=genotype(lMM, lMm, lmm, lN);
		igl++;
	} else {
		likelihoods.push_back(genotype(lMM, lMm, lmm, lN));	
	};
}

void population_genotypes::clear(){
	frozen=true;
	igl=likelihoods.begin();
}

genotype_pair convert(const population_genotypes& popgl, const count_t &a, count_t &b, const float_t &ma, const float_t &Ma, const float_t &mb, const float_t &Mb){
	if ((popgl.likelihoods[a].N>=ma) && (popgl.likelihoods[a].N<=Ma) && (popgl.likelihoods[b].N>=mb) && (popgl.likelihoods[b].N<=Mb) ){
		return genotype_pair (popgl.likelihoods[a].lMM, popgl.likelihoods[a].lMm, popgl.likelihoods[a].lmm, popgl.likelihoods[b].lMM, popgl.likelihoods[b].lMm, popgl.likelihoods[b].lmm, popgl.m);
	}
	else return genotype_pair (0, 0, 0, 0, 0, 0, 0);
}

std::string population_genotypes::header(void) const{
	std::string line="@ID0\tID1";
	std::vector <std::string>::const_iterator s_it=sample_names_.cbegin(), end=sample_names_.cend();
	while(s_it!=end){
		line+='\t';	
		line+=*s_it;
		s_it++;	
	}
	line+='\n';
	return line;
}

std::ostream& operator<< (std::ostream& out, const population_genotypes& x)
{
	out << x.id0 << '\t' << x.id1;
	std::vector <genotype>::const_iterator s_it=x.likelihoods.cbegin(), end=x.likelihoods.cend();
//	std::vector <genotype>                      likelihoods;		//!< genotypic likelihood
	while(s_it!=end){
		out << '\t' << *s_it;
		s_it++;	
	}
	return out;
}

std::istream& operator>> (std::istream& in, population_genotypes& x)
{
	return in;
}

std::ostream& operator<< (std::ostream& out, const genotype& x)
{
	out << x.lMM << '/' << x.lMm << '/' << x.lmm << '/' << x.N;
	return out;
}

std::istream& operator>> (std::istream& in, genotype& x)
{
	in >> x.lMM >> x.lMm >> x.lmm >> x.N;
	return in;
}
