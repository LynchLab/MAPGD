#include "data.h"

const std::string Data::file_name=".txt";//!< The destination table in the Db.
const std::string Data::table_name="NONE";//!< The destination table in the Db.
const bool Data::binary=false; //!< A flag to indicate that binary reading/writing is not supported by default.

static int nifty_counter=0; 

std::map <std::string, Data*(*)(const std::vector <std::string> &) > m_data_ctor;
/*static typename std::aligned_storage<sizeof (std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >), alignof (std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >)>::type
  m_data_ctor_buf; // memory for the stream object
std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >& stream = reinterpret_cast<std::map <std::string, Data*(*)(const std::vector <std::string>) >& > (m_data_ctor_buf);*/

/*std::vector <std::string> registry_list(void)
{
	std::vector<std::string> ret;
	for( std::map <std::string, Data*(*)(const std::vector <std::string> &)>::iterator it = m_data_ctor.begin(); it != m_data_ctor.end(); ++it) ret.push_back(it->first);
	return ret;
}*/

Registry_initalizer::Registry_initalizer()
{
	if (nifty_counter++ == 0) 
	{
#ifdef DEBUG
		std::cerr << "creating m_data_ctor\n";
#endif
		new (&m_data_ctor) std::map <std::string, Data*(*)(const std::vector <std::string> &)>; 
	}
}

Registry_initalizer::~Registry_initalizer()
{
	if(nifty_counter-- == 0) free(&m_data_ctor);
}

Data* Data::new_from_str(const std::string &Name, const std::vector<std::string> &columns)
{
	if ( m_data_ctor.find(Name) == m_data_ctor.end() ) {
		fprintf(stderr, gettext("mapgd:%s:%d: Cannot find class with name %s. Class names are:\n"), __FILE__, __LINE__, Name.c_str());
		for(std::map<std::string, Data*(*)(const std::vector <std::string> &) >::iterator it = m_data_ctor.begin(); it != m_data_ctor.end(); ++it) 
		{
  			std::cerr << '\t' << it->first << "\n";
		}	
		fprintf(stderr, gettext("mapgd:%s:%d: Exiting.\n"), __FILE__, __LINE__);
		exit(NOCLASS);
	} else 	{
		return m_data_ctor[Name](columns);
	}
}

void 
Data::sql_read(std::istream &in)
{
	std::string line;
	std::getline(in, line);
	std::cerr << line << std::endl;
}

std::istream& operator >> (std::istream& in, Data &data)
{
	data.read(in);
	return in;
}

std::ostream& operator << (std::ostream& out, const Data &data)
{
	data.write(out);
	return out;
}

Registration::Registration (const std::string &str, Data*(*fn)(const std::vector <std::string> &) )
{
#ifdef DEBUG
	std::cerr << "Registering " << str << std::endl;
	for(std::map<std::string, Data*(*)(const std::vector <std::string> &) >::iterator it = m_data_ctor.begin(); it != m_data_ctor.end(); ++it) 
	{
  		std::cerr << '\t' << it->first << "\n";
	}	
#endif
	m_data_ctor[str]=fn;
	std::vector <std::string> test;
#ifdef DEBUG
//	try 
//	{
		Data *test_data=m_data_ctor[str](test);
/*	}
	catch (const std::exception& e)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Class %s constructor throws exception %s. Exiting.\n"), __FILE__, __LINE__, str.c_str(), e.what().c_str() );
		exit(NOCLASS);
	}*/
#else
	Data *test_data=m_data_ctor[str](test);
#endif 

#ifdef DEBUG
	std::cerr << &m_data_ctor << ", " << test_data->get_table_name() << std::endl;
	for(std::map<std::string, Data*(*)(const std::vector <std::string> &) >::iterator it = m_data_ctor.begin(); it != m_data_ctor.end(); ++it) 
	{
  		std::cerr << '\t' << it->first << "\n";
	}	
#endif
	if (test_data->get_table_name()!=str)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Class %s does not return the correct type. Exiting.\n"), __FILE__, __LINE__, str.c_str() );
		exit(NOCLASS);
	}
	delete test_data;
	name_="FOOL!";//str;
}


Registration::~Registration ()
{
#ifdef DEBUG
	std::cerr << "Deregistering " << name_ << std::endl;
#endif
	m_data_ctor.erase(name_);
}


id1_t Indexed_data::get_abs_pos (void) const
{
	return abs_pos_;
}

id1_t Double_indexed_data::get_abs_pos (void) const
{
	return abs_pos1_;
}

id1_t Double_indexed_data::get_abs_pos2 (void) const
{
	return abs_pos2_;
}

void Indexed_data::set_abs_pos (const id1_t &new_pos)
{
	abs_pos_=new_pos;
}

void Double_indexed_data::set_abs_pos (const id1_t &new_pos)
{
	abs_pos1_=new_pos;
}

void Double_indexed_data::set_abs_pos2 (const id1_t &new_pos)
{
	abs_pos2_=new_pos;
}
