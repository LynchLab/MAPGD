#include "pooled_data.h"

const std::string Pooled_data::file_name=".pol";
const std::string Pooled_data::table_name="SAMPLE";
const Registration Pooled_data::registered=Registration(Pooled_data::table_name, Pooled_data::create);

Pooled_data::Pooled_data ()
{
	names_.clear();
	p.clear();
	polyll.clear();
	fixedll.clear();
	major.base=4;
	minor.base=4;
	abs_pos_=0;
	coverage=0;
	error=0;
	delim='\t';
}

void 
Pooled_data::set_sample_names (const std::vector <std::string> &columns)
{
        size_t size=columns.size();
        names_=columns;
        p.assign(size, 0);
        cov.assign(size, 0);
        polyll.assign(size, 0);
        fixedll.assign(size, 0);
}

Pooled_data::Pooled_data (const std::vector <std::string> &columns)
{
	size_t size=(columns.size()-6);
	names_=std::vector <std::string> (columns.cbegin()+6, columns.cend() );
//	names_.assign(size,"none");//=std::vector <std::string> (columns.cbegin()+6, columns.cend() );

        p.assign(size, 0);
        cov.assign(size, 0);
        polyll.assign(size, 0);
        fixedll.assign(size, 0);

        major.base=4;
        minor.base=4;
	abs_pos_=0;
        coverage=0;
        error=0;
        delim='\t';
}

void
Pooled_data::read (std::istream& in)
{
	std::string line, f;
	std::getline(in, line);
	std::stringstream line_stream(line);

	line_stream >> major;
	line_stream >> minor;
	line_stream >> coverage;
	line_stream >> error;
	line_stream.get();

	for (size_t s=0; s<names_.size();++s) {
		getline(line_stream, f, '/');
//		std::cerr << f << std::endl;
		if(f!="..."){
	        	p[s]=std::stof(f);
			getline(line_stream, f, '/');
			cov[s]=std::stof(f);
			getline(line_stream, f, '/');
			polyll[s]=std::stof(f);
			getline(line_stream, f, '\t');
			fixedll[s]=std::stof(f);
		} else {
			p[s]=NAN;
			getline(line_stream, f, '/');
			cov[s]=0;
			getline(line_stream, f, '/');
			polyll[s]=0;
			getline(line_stream, f, '\t');
			fixedll[s]=0;
		}
	}	
}

void
Pooled_data::write (std::ostream& out) const
{
	out << major << delim;
	out << minor << delim;
	out << coverage << delim;
	if (coverage!=0) out << error;
	else out << '.';

	for (size_t s=0; s<names_.size(); ++s) {
		if (!isnan(p[s]) ){
	//		out << std::setprecision(1);
			out << '\t' << p[s] << '/';
			out << cov[s] << '/';
			out << polyll[s] << '/';
			out << fixedll[s];
//			out << majorll[s];
		} else {
			out << "\t.../.../.../..";
		}
	}
}

std::string 
Pooled_data::header(void) const 
{
	std::string line="@SCFNAME       \tPOS\tMAJOR\tMINOR\tCOVRAG\tERROR";
	for (size_t s=0; s<names_.size(); ++s) {
		line+=('\t'+names_[s]);
	}
	line+='\n';
	return line;
}

size_t 
Pooled_data::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+names_.size()*sizeof(char);
}
	
Allele 
Pooled_data::to_allele(const size_t & x)
{
	Allele a;
        a.set_abs_pos(abs_pos_);             //!< the scaffold identifer of the allele.
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
	
Pooled_data & 
Pooled_data::operator=(const Pooled_data &rhs){
        this->major=rhs.major;
        this->minor=rhs.minor;
        this->coverage=rhs.coverage;
        this->error=rhs.error;
	
	this->names_=rhs.names_;
	this->p=rhs.p;
	this->cov=rhs.cov;
	this->polyll=rhs.polyll;
	this->fixedll=rhs.fixedll;

	return *this;	
}	
