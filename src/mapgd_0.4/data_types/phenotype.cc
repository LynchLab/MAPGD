#include "phenotype.h"

const std::string Phenotype::file_name=".phe";
const std::string Phenotype::table_name="PHENOTYPE";
const bool Phenotype::binary=false;

const Registration Phenotype::registered=Registration(Phenotype::table_name, Phenotype::create);

Phenotype::Phenotype ()
{
	delim='\t';
	n_traits_=0;
	n_samples_=0;
}

Phenotype::Phenotype (const std::vector<std::string> &fields)
{
	delim='\t';
	n_traits_=fields.size();
	n_samples_=0;
	if (n_traits_ > 1)
	{
		n_traits_--;
		trait=std::vector <std::string> (fields.begin()+1, fields.end() );
#ifdef EIGEN
//	value=Eigen::MatrixXf::Zero(0,n_traits_);
#else
	value=std::vector<std::vector <real_t> > (n_traits_);
#endif 
	}
}

Phenotype::Phenotype (const size_t &N)
{
	delim='\t';
	n_traits_=N;
	n_samples_=0;
#ifdef EIGEN
//	value=Eigen::MatrixXf::Zero(0,n_traits_);
#else
	value=std::vector<std::vector <real_t> > (n_traits_);
#endif 
	for (size_t x=0; x<N; x++)
	{
		trait.push_back(std::to_string(x) );
	}
}

void
Phenotype::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::vector <std::string> fields=split(line, delim);

	sample_name.push_back(fields[0]);
#ifdef EIGEN
	//value.conservativeResize(n_samples_++, n_traits_);
#endif 

	for (size_t y=0; y<n_traits_; y++)
	{
#ifdef EIGEN
		//value(n_samples_-1,y)=atof(fields[1+y].c_str());
#else
		value[y].push_back(atof(fields[1+y].c_str()));
#endif 
	}
}

void
Phenotype::write (std::ostream& out) const 
{
	for (size_t x=0; x<n_samples_; x++)
	{
		std::cout << sample_name[x] << delim;
		for (size_t y=0; y<n_traits_; y++)
		{
#ifdef EIGEN
	//		std::cout << value(x,y) << delim;
#else
			std::cout << value[y][x] << delim;
#endif
		}
		if (x!=n_samples_-1) std::cout << std::endl;
	}
}

void
Phenotype::add_sample(const uint32_t &id, const real_t *new_value)
{
	sample_name.push_back( std::to_string(id) );
	for (size_t y=0; y<n_traits_; y++)
	{
#ifdef EIGEN
//	std::cout << value(x,y) << delim;
#else
	value[y].push_back(new_value[y]);
#endif
	}
	n_samples_++;
}

std::string 
Phenotype::header(void) const 
{
	std::stringstream ret; 
	ret << "@SMPNAME:" << n_samples_;
	for (size_t x=0; x<n_traits_; x++)
	{
		 ret << delim << trait[x];
	}
	ret << std::endl;
	return ret.str();
}


size_t 
Phenotype::size(void) const 
{
	//TODO IMPLEMENT
	return 0;
}

/*
void 
Phenotype::clear(void)
{
	zero();
}

void 
Phenotype::zero(void)
{
	//TODO IMPLEMENT
}*/

/*
const size_t
Phenotype::sample_size(void) const
{
	return n_samples_;
}*/;

const bool 
Phenotype::get_binary(void) const
{
	return binary;
}

const std::string Phenotype::get_file_name(void) const
{
	return file_name;
}

const std::string Phenotype::get_table_name(void) const
{
	return table_name;
}
