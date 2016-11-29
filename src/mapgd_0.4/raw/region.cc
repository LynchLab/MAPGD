#include "region.h"
Region::Region()
{
	scf_name="";
	start=0;
	stop=0;
}

Region::Region(const std::string &this_name, const id1_t &this_start, const id1_t &this_stop)
{
	scf_name=this_name;
	start=this_start;
	stop=this_stop;
}

Region::Region(const id1_t &this_start, const id1_t &this_stop)
{
	scf_name="";
	abs_start=this_start;
	abs_stop=this_stop;
}

Region::Region(const id0_t &this_id0, const id1_t &this_start, const id1_t &this_stop)
{
	id0=this_id0;
	start=this_start;
	stop=this_stop;
}

void
Region::set(const File_index &index)
{
	if (scf_name!="") 
	{
		abs_start=index.get_abs_pos(scf_name, start);
		abs_stop=index.get_abs_pos(scf_name, stop);
	} else {
		abs_start=0;
		abs_stop=index.get_reference_size()+1;
	}
}

std::ostream& operator << (std::ostream& out, const Region& x)
{
	out << x.scf_name << ":" << x.start << "-" << x.stop;
	return out;
}

std::istream& operator >> (std::istream& in, Region& x)
{
	std::string line;
	getline(in, line, ':');
	x.scf_name=line;
	getline(in, line, '-');
	x.start=std::stoi(line);
	getline(in, line, '\t');
	x.stop=std::stoi(line);
	return in;
}

Region ator(const char *c_str)
{
	Region region;
	std::string str(c_str);
	std::vector <std::string> feilds=split(str, ':');
	if (feilds.size()==2)
	{
		region.scf_name=feilds[0];
		feilds=split(feilds[1],'-');
		region.start=std::stoi(feilds[0]);
		region.stop=std::stoi(feilds[1]);
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << "error parsing region.\n";
	}
	return region;
}

bool isregion(const char *c_str)
{
	return true;	
}
