#include "locus.h"

const std::string Locus::file_name=".pro";
const std::string Locus::table_name="QUARTETS";
const bool Locus::binary=true;

const Registration Locus::registered=Registration(Locus::table_name, Locus::create);
const gt_t Locus::default_order[5] = {0,1,2,3,4};
	
/*! \brief Locus is initialized to abs_pos_=0, id1=0, and set to contain 0 samples.
*/	
Locus::Locus(void)
{
	sample.clear();
	sample_names_.clear();
	abs_pos_=0;
	ref.base=4;
	all_unredacted_=true;
	redacted_samples_.clear();
}

/*! \brief Locus is initialized to abs_pos_=0, id1=0, and set to contain 0 samples.
*/	
Locus::Locus(const count_t &size)
{
	sample.assign(size, quartet() );
	sample_names_.clear();//=std::vector<std::string> (size);
	char buffer[21];
	for (size_t x=0; x<size; ++x){
		snprintf(buffer, 20, "Sample_%d", x+1);
		sample_names_.push_back(std::string(buffer) );
	}

	all_unredacted_=true;
	redacted_samples_.resize(size);
	std::fill(redacted_samples_.begin(), redacted_samples_.end(),false);
	
	abs_pos_=0;
	ref.base=4;
}

Locus::Locus(const std::vector<std::string> &column_names)
{
	if (column_names.size()>=4)
	{
		//
		sample_names_=std::vector<std::string> (&column_names[3], &column_names.back()+1 );
		sample.assign(sample_names_.size(), quartet() );
	}
	else 
	{
		sample.clear();
		sample_names_.clear();
	}
	abs_pos_=0;
	ref.base=4;

	all_unredacted_=true;
	redacted_samples_.resize(sample_names_.size() );
	std::fill(redacted_samples_.begin(), redacted_samples_.end(),false);
}

Locus & 
Locus::operator =(const Locus& rhs)
{
	memcpy(sorted_, Locus::default_order, 5*sizeof(gt_t) );
	sample.clear();
        sample=rhs.sample; 
	sample_names_.clear();
	sample_names_=rhs.sample_names_;
	redacted_samples_.clear();
	redacted_samples_=rhs.redacted_samples_;
	all_unredacted_=rhs.all_unredacted_;
        abs_pos_=rhs.abs_pos_;
	return *this;
}

Locus & 
Locus::operator+=(const Locus &rhs) 
{
	if (rhs.abs_pos_==this->abs_pos_ ){
		sample.reserve( sample.size() + rhs.sample.size() ); // preallocate memory
		sample.insert( sample.end(), rhs.sample.begin(), rhs.sample.end() );
		sample_names_.reserve(sample_names_.size() + rhs.sample_names_.size() );
		sample_names_.insert( sample_names_.end(), rhs.sample_names_.begin(), rhs.sample_names_.end() );
		redacted_samples_.reserve(redacted_samples_.size()+rhs.redacted_samples_.size() );
		redacted_samples_.insert(redacted_samples_.end(), rhs.redacted_samples_.begin(), rhs.redacted_samples_.end() );
		return *this;
	}
	std::cerr << __FILE__ << ":" << __LINE__ << "operator undefined for different loci.\n";
	return *this;
}

const 
Locus Locus::operator +(const Locus& rhs) const 
{
	return Locus (*this) += rhs;
}

/*
void Locus::unmask(quartet_t *q){
	q->masked=false;
}

void Locus::mask(quartet_t *q){
	q->masked=true;
}

void Locus::unmask(count_t a){
	if(a<sample.size() ) sample[a].masked=false;
}*/


/// Unmask all the quartets at this locus.
void 
Locus::unmaskall(void)
{
	for (size_t s=0; s<sample.size();++s){
		sample[s].masked=false;
	}
}


/// Mask all the quartets at this locus.
void 
Locus::maskall(void)
{
	for (size_t s=0; s<sample.size();++s){
		sample[s].masked=true;
	};
}

/// Unmask all the quartets in the populations in the vector.
void 
Locus::unmask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<s.size(); ++s_){
		if (s_>=sample.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to unmask a non-existent quartet. Exiting."; exit(0); };
		sample[ s[s_] ].masked=false;
	}
}


/// Mask all the quartets in the populations in the vector.
void 
Locus::mask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<s.size();++s_){
		if (s_>=sample.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to mask a non-existent quartet. Exiting."; exit(0); };
		sample[s_].masked=true;
	};
}

count_t 
Locus::maskedcount(void) const
{
	count_t count=0;
	for(std::vector<quartet>::const_iterator it = sample.begin(); it != sample.end(); ++it) if (it->masked) count++;
	return count;
}

void 
Locus::sort(void)
{
	count_t total[5]={0};
	
	for (unsigned int s=0; s<sample.size();++s){
		if (sample[s].masked) continue;
		total[0]+=sample[s].base[0];
		total[1]+=sample[s].base[1];
		total[2]+=sample[s].base[2];
		total[3]+=sample[s].base[3];
		total[4]+=sample[s].base[4];
	};
	if (total[sorted_[0]]<total[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (total[sorted_[1]]<total[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (total[sorted_[2]]<total[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (total[sorted_[0]]<total[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (total[sorted_[1]]<total[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
}

void 
Locus::sort(count_t s)
{
	if (sample[s].base[sorted_[0]]<sample[s].base[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (sample[s].base[sorted_[1]]<sample[s].base[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (sample[s].base[sorted_[2]]<sample[s].base[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (sample[s].base[sorted_[0]]<sample[s].base[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (sample[s].base[sorted_[1]]<sample[s].base[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
}

void 
Locus::swap(count_t x, count_t y)
{	
	std::swap(sorted_[x], sorted_[y]);
}

count_t 
Locus::getcount(count_t s, count_t c) const
{
	if (c<5) return sample[s].base[sorted_[c]];
	else return 0;
}

const count_t 
Locus::getindex(count_t c) const
{
	if (c<5) return  sorted_[c];
	else return -1;
}

count_t 
Locus::getcount(count_t c) const
{
	count_t total=0;
	if (c<5){
		for (unsigned int s=0; s<sample.size() ;++s){
			if (sample[s].masked) continue;
			total+=sample[s].base[sorted_[c]];
		}
		return total;
	}
	else return 0;
}

count_t 
Locus::getcoverage(count_t s) const
{
	if (s<sample.size() ) return sample[s].base[0]+
				sample[s].base[1]+
				sample[s].base[2]+
				sample[s].base[3];
	else {
		std::cerr << "mapgd:" << __FILE__ << ":" << __LINE__ << ": Attempted to access a non-existent sample." << std::endl;
		exit(0);
	};
}

count_t 
Locus::getcoverage() const
{
	count_t total=0;
	for (size_t s=0; s<sample.size();++s){
		if (sample[s].masked) continue;
		total+=sample[s].base[0]+
			sample[s].base[1]+
			sample[s].base[2]+
			sample[s].base[3];
	};
	return total;
}

count_t 
Locus::get_max_coverage() const
{
	count_t max_dp=0;
	for (size_t s=0; s<sample.size();++s){
		if (sample[s].masked) continue;
		max_dp = std::max(count_t(sample[s].base[0]+
			sample[s].base[1]+
			sample[s].base[2]+
			sample[s].base[3]), max_dp);
	};
	return max_dp;
}

count_t 
Locus::getcoverage_nomask() const
{
	count_t total=0;
	for (size_t s=0; s<sample.size();++s){
		total+=sample[s].base[0]+
			sample[s].base[1]+
			sample[s].base[2]+
			sample[s].base[3];
	};
	return total;
}

/// sets quartet_t c to q.
void 
Locus::set_quartet(const quartet_t &q, const count_t &c)
{
	sample[c]=q;
}		

/// returns quartet_t c.
const quartet_t & 
Locus::get_quartet(const count_t &c) const 
{
	return sample[c];
}

/// returns quartet_t c.
quartet_t & 
Locus::get_quartet(const count_t &c) 
{
	return sample[c];
}
	
/// sets the number of sampels to c.
void 
Locus::resize(const size_t &c)
{
	sample.resize(c);
	sample_names_.resize(c);
	redacted_samples_.resize(c);
}

void
Locus::mask_low_cov( const count_t &dp )
{
	//The parens are unneccisary, but I don't like getting compiler warnings, so.... 
	for (size_t s=0; s<sample.size();++s) sample[s].masked=(sample[s].masked| ( count(sample[s])<=dp ) );
}

void
Locus::write (std::ostream& out) const
{
	std::vector <bool>::const_iterator r=redacted_samples_.cbegin();
	out << ref;
	for (std::vector <quartet_t>::const_iterator s=sample.cbegin(); s<sample.cend(); s++) {
		if (!*(r++) ) 
			out << '\t' << *s;
	}
}

void
Locus::write_binary (std::ostream& out) const
{
	out.write((char *)&ref.base, sizeof(gt_t) );
	if (all_unredacted_)
		out.write((char *)&sample[0], (size_t)(sample.size()*sizeof(quartet_t) ) );
	else {
		for (std::vector <quartet_t>::const_iterator s=sample.cbegin(); s<sample.cend(); s++) 
			out.write( (char *)&(*s), (size_t)(sizeof(quartet_t) ) );
	}
}

void
Locus::read_binary (std::istream& in) 
{
	in.read((char *)&ref.base, sizeof(gt_t) );
	in.read((char *)&sample[0], (size_t)(sample.size()*sizeof(quartet_t) ) );
}

void
Locus::read (std::istream &in)
{
	//Set sorted_ back to default (a,c,g,t)
	memcpy(sorted_, Locus::default_order, 5*sizeof(gt_t) );

	std::string line;
	//remove characters from the in stream until a new line.
	std::getline(in, line);
	std::stringstream line_stream(line);
	
	line_stream >> ref;

	for (std::vector <quartet_t>::iterator s=sample.begin(); s<sample.end(); s++) {
		line_stream >> *s;
	}
	//We do not save mask states, so all data read in is unmasked
	unmaskall();				
}

std::istream &
mpileup (std::istream& in, Locus& x, const int &offset, const int &ncolumns)
{
	std::string line;
	std::vector <std::string> column;
	std::vector <std::string>::iterator column_it=column.begin();

	in.get();
	getline(in, line);
	column=split(line, '\t');

	if (column.size()!=x.sample.size()*ncolumns+offset-1) 
	{
		std::cerr << __FILE__ << "." << __LINE__ << ": unexpected error parsing mpileup file." << std::endl;
		std::cerr << column.size() << " not " << x.sample.size()*ncolumns+offset-1 << std::endl;
		in.clear (in.rdstate() | std::ios::failbit);
		return in;
	}

	x.ref=Base::ctob(column[0].c_str()[0] );
	
	for (size_t s=0; s<x.sample.size(); ++s) {
		memset(x.sample[s].base, 0, sizeof(count_t)*5 );
		scan(x, column[s*ncolumns+offset], x.sample[s] );
	}
	return in;
}
	
std::string 
Locus::header(void) const
{

	std::string line="@ID0       \tID1\tREF";
	for (size_t s=0; s<sample_names_.size(); ++s) {
		if (!redacted_samples_[s] ) 
			line+=('\t'+sample_names_[s]);
	}
	line+='\n';
	return line;
}

size_t
Locus::size(void) const
{
	return sizeof(count_t)*(5)+sizeof(quartet_t)*sample.size()+sizeof(id1_t);
}

void scan(const Locus & site, const std::string &str, quartet_t &q)
{
	size_t skip;
	std::string number;
	std::string::const_iterator it=str.begin();
	std::string::const_iterator end=str.end();
	size_t lim=0;
	while (it!=end) {
		if (*it == '^') {it++; if(it==end) break; it++; continue;}
		if (*it == '$' || *it=='*') {it++; continue;}
		else if (*it == '+' || *it == '-') {
			number = "";
			it++;
			while(isdigit(*it) ) {
				number.push_back(*it); 
				it++; 
			}
			lim=atoi(number.c_str() );
 			for(skip=0; skip<lim; ++skip) it++;
			continue;
    		} else if(*it == ',' || *it=='.') {
			q.base[site.ref.base]++;
			#ifdef DEBUG
				if (site.ref.base==5) {
					std::cerr << __FILE__ << ":" << __LINE__ << ": no reference nucleotide set.\n";
				}
			#endif 
		} else {
			q.base[Base::ctob(*it)]++;
		}
		it++;
	}
}

const std::string 
Locus::sql_header(void) const {
	return "(ABS_POS int, SMPNUM int, A int, C int, G int, T int, PRIMARY KEY (ABS_POS, SMPNUM) )";
}

const std::string 
Locus::sql_column_names(void) const {
	return "(ABS_POS, SMPNUM, A, C, G, T)";
}

const std::string 
Locus::sql_get_constructor(void) const {
	return "SELECT SMPNAME FROM SAMPLE";
}

void
Locus::sql_read(std::istream &in)
{
	//id1 last_pos=-1;
	count_t smpnum;
	size_t size=sample.size();
	for (size_t x=0; x<size; ++x)
	{
		in >> abs_pos_;
		in >> smpnum;
		in >> sample[x].base[0];
		in >> sample[x].base[1];
		in >> sample[x].base[2];
		in >> sample[x].base[3];
	} 
	unmaskall();				
} 

const std::string 
Locus::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
	char *write_ptr=return_buffer;
	std::vector<quartet>::const_iterator it=cbegin(), end=cend();
	int x=0;

	write_ptr+=snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %d, %d, %d, %d, %d)", 
       		abs_pos_,
		++x,
		(*(it++))[0],
		(*it)[1],
		(*it)[2],
		(*it)[3]);
	while (it!=end) {
		write_ptr+=snprintf(write_ptr, SQL_LINE_SIZE, ", (%d, %d, %d, %d, %d, %d)", 
       		abs_pos_,
		++x,
		(*(it++))[0],
		(*it)[1],
		(*it)[2],
		(*it)[3]);
        }
	return std::string(return_buffer);
}

const std::string
Locus::get_file_name(void) const
{
	return Locus::file_name;
}

const std::string
Locus::get_table_name(void) const
{
	return Locus::table_name;
}

const bool 
Locus::get_binary(void) const
{
	return binary;
}

void 
Locus::unredact_all(void)
{
	all_unredacted_=true;
	for (int x=0;x<redacted_samples_.size();++x) {
		redacted_samples_[x]=false;
	}
}

void 
Locus::redact_all(void)
{
	all_unredacted_=false;
	for (int x=0;x<redacted_samples_.size();++x) {
		redacted_samples_[x]=true;
	}
}

void 
Locus::unredact(const size_t &x)
{
	//if ( ) {
		redacted_samples_[x]=false;
//		all_unredacted_=all_unredacted();
	/*}
	else {
	}*/
}

void 
Locus::redact(const size_t &x)
{
	all_unredacted_=false;
	//if ( )
		redacted_samples_[x]=true;
	//else {
	//}	
}
