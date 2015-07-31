#include "locus.h"

locus::locus(void){};

///This form of the constructor initalizes locus with size samples.
locus::locus(size_t size){
	sample.assign(size, quartet() );
	sorted_=memcopy(?,?);		
}

locus & locus::operator =(const locus& rhs)
{
        sample=rhs.sample; 	
	sorted_=memcopy(?,?);		
	return *this;
}

///NOTE: the order of the locus is preserved through this operation.
locus & locus::operator+=(const locus &rhs) 
{
	sample.reserve( sample.size() + rhs.sample.size() ); // preallocate memory
	sample.insert( sample.end(), rhs.sample.begin(), rhs.sample.end() );
	return *this;
}

///NOTE: the order of the locus is preserved through this operation.
const locus locus::operator +(const locus& rhs) const {
	return locus (*this) += rhs;
}

/// Unmask all the quartets at this locus.
void locus::unmaskall(void)
{
	for (size_t s=0; s<sample.size();++s){
		sample[s].masked=false;
	}
}


/// Mask all the quartets at this locus.
void locus::maskall(void)
{
	for (size_t s=0; s<sample.size();++s){
		sample[s].masked=true;
	};
}

/// Unmask all the quartets in the populations in the vector.
void locus::unmask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<sample.size();++s_){
		if (s_>=sample.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to unmask a non-existent quartet. Exiting."; exit(0); };
		sample[s_].masked=false;
	}
}


/// Mask all the quartets in the populations in the vector.
void locus::mask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<s.size();++s_){
		if (s_>=sample.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to mask a non-existent quartet. Exiting."; exit(0); };
		sample[s_].masked=true;
	};
}

size_t locus::maskedcount(void) const
{
	size_t count=0;
	for(std::vector<quartet>::const_iterator it = sample.begin(); it != sample.end(); ++it) if (it->masked) count++;
	return count;
}

void locus::sort(void)
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

void locus::sort(size_t s)
{
	if (s<sample.size() ) {
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
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted to access a non-existent sample." << std::endl;
		exit(0);
	}
}

void locus::swap(const gt_t &lhs, const gt_t &rhs)
{	
	std::swap(sorted_[lhs], sorted_[rhs]);
}

count_t locus::get_count(count_t s, count_t c) const
{
	if (c<5) return sample[s].base[sorted_[c]];
	else return 0;
}

const count_t locus::get_index(count_t c) const
{
	if (c<5) return  sorted_[c];
	else return -1;
}

count_t locus::get_count(count_t c) const
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

count_t locus::get_coverage(count_t s) const
{
	if (s<sample.size() ) return sample[s].base[0]+
				sample[s].base[1]+
				sample[s].base[2]+
				sample[s].base[3];
	else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted to access a non-existent sample." << std::endl;
		exit(0);
	};
}

count_t locus::get_coverage() const
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

/// sets quartet_t c to q.
void locus::set_quartet(const quartet_t &q, const count_t &c)
{
	sample[c]=q;
}		

/// returns quartet_t c.
const quartet_t & locus::get_quartet(const count_t &c) const 
{
	return sample[c];
}

/// returns quartet_t c.
quartet_t & locus::get_quartet(const count_t &c) 
{
	return sample[c];
}
	
/// no return type.
void locus::resize(const size_t &c)
{
	sample.resize(c);
}

/// no return type. TODO User iterator?
void locus::mask_low_cov( const count_t &dp )
{ 
	for (size_t s=0; s<sample.size();++s) sample[s].masked=(count(sample[s])<dp);
}
