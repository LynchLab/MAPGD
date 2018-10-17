#ifndef _PLINK_PHENO_H_
#define _PLINK_PHENO_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "external-file.h"
#include "external-data.h"

#include "typedef.h"
#include "datatypes.h"
#include "phenotype.h"
#include "error_codes.h"

enum Binary { missing, unaffected, affected};
enum Sex { unknown, male, female};

template <typename T> 
class Plink_record  
{
	public:
	uint32_t family_id;
	uint32_t individual_id;
	uint32_t maternal_id;
	uint32_t paternal_id;
	size_t n_values;
	Sex sex;
	T *values;

	Plink_record()
	{
		n_values=0;
		values=NULL;
	}

	Plink_record(size_t n)
	{
		n_values=n;
		values=new T[n_values];
	}

	void
	resize(size_t n)
	{
		if(n_values!=n)
		{
			if (values!=NULL)
			{
				delete values;
				values=NULL;
			}
			n_values=n;
			values=new T[n_values];
		}
	}

	~Plink_record()
	{
		if (values!=NULL)
		{
			delete values;
			values=NULL;
		}
	}

	Plink_record& operator =(const Plink_record& rhs)
	{
		if (&rhs != this) 
		{
			family_id=rhs.family_id;
			individual_id=rhs.individual_id;
			maternal_id=rhs.maternal_id;
			paternal_id=rhs.paternal_id;
			sex=rhs.sex;
			if (n_values!=rhs.n_values) 
			{
				if (values!=NULL)
				{
					delete values;
				}
				values=new T[n_values];
				n_values=rhs.n_values;
			}
			memcpy(values, rhs.values, n_values*sizeof(T) );
		}
		return *this;
	}
};

class Plink_data : public External_data {
private:
	Plink_record <real_t> record_;
	std::vector <std::string > sample_names_;
	real_t missing_value;
	size_t n_samples_;
public:
	Plink_data(const size_t &);

	Plink_data(){};

	void set_family_id(const std::string &);
	void set_individual_id(const std::string &);
	int set_values(const std::vector <std::string> &);
	int set_values(const real_t *);
	size_t trait_size(void) const;

	void put (const Data *, ...);
	void get (Data *, ...) const;

	void put (const Phenotype &);
	void get (Phenotype &) const;

	std::vector<std::string> get_sample_names (void) const;
};

//'FID' and 'IID

class Plink_file : public External_file <Plink_data> {
private:
	using Base_file::in_;
	using Base_file::out_;
	using Base_file::open_;
	using Base_file::open_no_extention;
	using Base_file::delim_column_;
public:
	void open(const std::ios_base::openmode &);	
	void open(const char *, const std::ios_base::openmode &);
	void close(void);
	void read (Plink_data &);
	void write (const Plink_data &);
	void write_header (const Plink_data &);
	Plink_data read_header();
};

#endif
