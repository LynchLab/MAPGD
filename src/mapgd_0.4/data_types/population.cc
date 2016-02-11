#include "population.h"

const std::string Population::file_name="gcf";
const std::string Population::table_name="GENOTYPES";
const Registration Population::registered=Registration(Population::table_name, Population::create);
/** @breif constuctor w/ initial values. **/

Population::Population(const std::vector <std::string> &column_names)
{
	sample_names_=std::vector <std::string> (column_names.cbegin()+3, column_names.cend() );
	likelihoods.resize(sample_names_.size() );
}

Population::Population()
{
}

Population::Population(const Population &rhs)
{
	*this=rhs;
}

Population::~Population(){}

/**@breif return size of Population if Population is set, 0 otherwise**/
size_t 
Population::size() const
{
	return likelihoods.size();
}

Population & 
Population::operator= (const Population& rhs)
{
	sample_names_=rhs.sample_names_;	//!< a vector of sample names.
	likelihoods=rhs.likelihoods;		//!< genotypic likelihood
	igl_=likelihoods.begin();		//!< an iterator to allow us to iterate over the likelihoods.

	major=rhs.major;			//!< identity of the major allele
	minor=rhs.minor;			//!< identity of the minor allele
	m=rhs.m;				//!< minor allele frequency
	abs_pos_=rhs.abs_pos_;			//!< scaffold number
	return *this;
}

std::string 
Population::header(void) const
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

void
Population::write (std::ostream& out) const
{
	std::vector <Genotype>::const_iterator s_it=likelihoods.cbegin(), end=likelihoods.cend();
	while(s_it!=end){
		out << '\t' << *s_it;
		s_it++;	
	}
}

void
Population::read (std::istream& in)
{
        std::string line;
        std::getline(in, line);
        std::stringstream line_stream(line);

	std::vector <Genotype>::iterator s_it=likelihoods.begin(), end=likelihoods.end();
	while(s_it!=end){
		line_stream >> *s_it;
		s_it++;
	}	
}
