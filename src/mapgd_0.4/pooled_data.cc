#include "pooled_data.h"

const std::string Pooled_data::file_name=".pol";
const std::string Pooled_data::table_name="SAMPLE";

Pooled_data::Pooled_data ()
{
	names_.clear();
	p.clear();
	polyll.clear();
	majorll.clear();
	minorll.clear();
	major.base=4;
	minor.base=4;
	id0=0;
	id1=0;
	coverage=0;
	error=0;
	delim='\t';
}

void Pooled_data::set_sample_names (const std::vector <std::string> &columns)
{
        size_t size=columns.size();
        names_=columns;

        p.assign(size, 0);
        polyll.assign(size, 0);
        majorll.assign(size, 0);
        minorll.assign(size, 0);
}

Pooled_data::Pooled_data (const std::vector <std::string> &columns)
{
	size_t size=(columns.size()-6);
	names_=std::vector<std::string> (size);

        p.assign(size, 0);
        polyll.assign(size, 0);
        majorll.assign(size, 0);
        minorll.assign(size, 0);

        major.base=4;
        minor.base=4;
        id0=0;
        id1=0;
        coverage=0;
        error=0;
        delim='\t';
}

std::istream& operator >> (std::istream& in, Pooled_data& x)
{
	std::string line;
	std::getline(in, line);
	std::stringstream line_stream(line);

	line_stream >> x.id1;
	line_stream >> x.major;
	line_stream >> x.minor;
	line_stream >> x.coverage;
	line_stream >> x.error;

	for (size_t s=0; s<x.names_.size();++s) {
		line_stream >> x.p[s];
		line_stream >> x.polyll[s];
		line_stream >> x.minorll[s];
		line_stream >> x.majorll[s];
	}
	return in;
}

std::ostream& operator<< (std::ostream& out, const Pooled_data& x) 
{
	out << x.id1 << x.delim;
	out << x.major << x.delim;
	out << x.minor << x.delim;
	out << x.coverage << x.delim;
	if (x.coverage!=0) out << x.error << x.delim;
	else out << '.' << x.delim;

	for (size_t s=0; s<x.names_.size(); ++s) {
		if (!isnan(x.p[s]) ){
	//		out << std::setprecision(1);
			out << x.p[s] << '/';
			out << x.polyll[s] << '/';
			out << x.minorll[s] << '/';
			out << x.majorll[s] << '\t';
		} else {
			out << ".../.../.../...\t";
		}
	}
	return out;
}

std::string Pooled_data::header(void) const 
{
	std::string line="@SCFNAME       \tPOS\tMAJOR\tMINOR\tCOVRAG\tERROR";
	for (size_t s=0; s<names_.size(); ++s) {
		line+=('\t'+names_[s]);
	}
	line+='\n';
	return line;
}

size_t Pooled_data::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+names_.size()*sizeof(char);
}
	
Allele Pooled_data::to_allele(const size_t & x)
{
	Allele a;
        a.id0=id0;             //!< the scaffold identifer of the allele.
        a.id1=id1;             //!< the bp location of the allele.
        a.freq=p[x];           //!< frequency of major allele.

        a.ref=4;               //!< idenity of ref allele.
        a.minor=minor.base;    //!< idenity of minor allele.
        a.major=major.base;    //!< idenity major allele.
        a.e1=4;                //!< idenity of error1
        a.e2=4;                //!< idenity of error2.

        a.error=error;         //!< ml error rate.
	a.coverage=coverage;   //!< population coverage.
	return a;	
}
