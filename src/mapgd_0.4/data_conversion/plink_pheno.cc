#include "plink_pheno.h"

Plink_data::Plink_data (const size_t &n)
{
	record_.resize(n);
}

void 
Plink_file::read (Plink_data &plink)
{
	if (open_) 
	{
		std::string line;
		std::getline(*in_, line);
		if (line.size()==0) 
		{
			close_table();
			return;
		}
		std::vector <std::string> fields=split(line, delim_column_);
		if (fields.size() > 2 )
		{
			plink.set_family_id( fields[0]);
			plink.set_individual_id( fields[1]);
			fields.erase(fields.begin(), fields.begin() + 2);
			if( plink.set_values( fields ) == 0) return;
			fprintf(stderr, gettext("mapgd:%s:%d: attempt to read from file failed. Expected %d columns, got %d.\n"), __FILE__, __LINE__, plink.trait_size()+2, fields.size()+2 );
			close_table();
		}
		fprintf(stderr, gettext("mapgd:%s:%d: attempt to read from file failed. Expected %d columns, got %d.\n"), __FILE__, __LINE__, plink.trait_size()+2, fields.size() );
		close_table();
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: attempt to read from unopened file.\n"), __FILE__, __LINE__ );
		exit(0);
	}
}

void
Plink_data::set_family_id(const std::string &str)
{
	record_.family_id=atoi(str.c_str() );
}

void
Plink_data::set_individual_id(const std::string &str)
{
	record_.individual_id=atoi(str.c_str() );
}

//TODO unwrap this or something?
int
Plink_data::set_values(const std::vector<std::string> &str)
{
	if(str.size()==record_.n_values)
	{
		for (size_t x=0; x<record_.n_values; x++)
		{
			record_.values[x]=atof(str[x].c_str() );
		}
		return 0;
	} else return BADSIZE;
}

void
Plink_data::get(Phenotype &p) const
{
	p.add_sample(record_.individual_id, record_.values);
}


void 
Plink_data::get(Data *data, ...) const
{
	if (0)
	{
	//	get_text("GOOD WORK %s", "BOO!");
	} else {
	//	get_text("BAD WORK %s", "BOO!");
	}
}

void
Plink_data::put(const Data *data, ...)
{
	if (0)
	{
	//	get_text("GOOD WORK %s", "BOO!");
	} else {
	//	get_text("BAD WORK %s", "BOO!");
	}
}

void
Plink_file::open(const std::ios_base::openmode &mode)
{
	fprintf(stderr, gettext("mapgd:%s:%d: cannot open file in std::ios_base::openmode &mode=%d.\n"), __FILE__, __LINE__, mode );
}

void
Plink_file::open(const char *file_name, const std::ios_base::openmode &mode)
{
	open_no_extention(file_name, mode);
/*
	if (mode & READ)
	{
		file_.open(file_name, READ);
		if ( file_) open_=true;
		else {
			fprintf(stderr, gettext("mapgd:%s:%d: Could not open file \"%s\" for reading.\n"), __FILE__, __LINE__, file_name );
			exit(EXIT_FAILURE);
		}
	}
	else if (mode & WRITE)
	{
		file_.open(file_name, WRITE);
		if ( file_) open_=true;
	}

	if (open_)
	{
		if (mode & READ)
		{
		}
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Unexpected error opening file \"%s\".\n"), __FILE__, __LINE__, file_name );
	}
*/		
}

void 
Plink_file::close(void)
{
	if (open_)
	{
		file_.close();
		open_=false;
	}
		
}

std::vector<std::string> 
Plink_data::get_sample_names (void) const 
{
	fprintf(stderr, gettext("mapgd:%s:%d: Retrieving sample names currently unsupported.\n"), __FILE__, __LINE__ );
}

size_t 
Plink_data::trait_size(void) const
{
	return record_.n_values;
}

Plink_data
Plink_file::read_header(void)
{
	int n=1;
	Plink_data plink(n);
	table_open_=true;
	return plink;
}
