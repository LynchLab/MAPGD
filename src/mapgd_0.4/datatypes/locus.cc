#include "locus.h"

locus::locus(void){
	locus(0);
};

///This form of the constructor initalizes locus with size sample_s.
locus::locus(const size_t &size){
	sample_.assign(size, quartet() );
	memcpy(sorted_, locus::defaultorder, sizeof(gt_t)*5);		
}

locus & locus::operator =(const locus& rhs)
{
        sample_=rhs.sample_; 	
	memcpy(sorted_, rhs.sorted_, sizeof(gt_t)*5);		
	return *this;
}

///NOTE: the order of the locus is preserved through this operation.
locus & locus::operator+=(const locus &rhs) 
{
	sample_.reserve( sample_.size() + rhs.sample_.size() ); // preallocate memory
	sample_.insert( sample_.end(), rhs.sample_.begin(), rhs.sample_.end() );
	return *this;
}

///NOTE: the order of the locus is preserved through this operation.
const locus locus::operator +(const locus& rhs) const {
	return locus (*this) += rhs;
}

/// Unmask all the quartets at this locus.
void locus::unmaskall(void)
{
	for (size_t s=0; s<sample_.size();++s){
		sample_[s].masked=false;
	}
}


/// Mask all the quartets at this locus.
void locus::maskall(void)
{
	for (size_t s=0; s<sample_.size();++s){
		sample_[s].masked=true;
	};
}

/// Unmask all the quartets in the populations in the vector.
void locus::unmask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<sample_.size();++s_){
		if (s_>=sample_.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to unmask a non-existent quartet. Exiting."; exit(0); };
		sample_[s_].masked=false;
	}
}


/// Mask all the quartets in the populations in the vector.
void locus::mask(const std::vector <size_t> &s)
{
	for (size_t s_=0; s_<s.size();++s_){
		if (s_>=sample_.size() ) {std::cerr << __FILE__ << ":" << __LINE__ << ":attempted to mask a non-existent quartet. Exiting."; exit(0); };
		sample_[s_].masked=true;
	};
}

size_t locus::maskedcount(void) const
{
	size_t count=0;
	for(std::vector<quartet>::const_iterator it = sample_.begin(); it != sample_.end(); ++it) if (it->masked) count++;
	return count;
}

void locus::sort(void)
{
	count_t total[5]={0};
	
	for (unsigned int s=0; s<sample_.size();++s){
		if (sample_[s].masked) continue;
		total[0]+=sample_[s].base[0];
		total[1]+=sample_[s].base[1];
		total[2]+=sample_[s].base[2];
		total[3]+=sample_[s].base[3];
		total[4]+=sample_[s].base[4];
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

void locus::sort(const size_t &s)
{
	if (s<sample_.size() ) {
		if (sample_[s].base[sorted_[0]]<sample_[s].base[sorted_[2]])
			std::swap(sorted_[0], sorted_[2]);
		if (sample_[s].base[sorted_[1]]<sample_[s].base[sorted_[3]])
			std::swap(sorted_[1], sorted_[3]);
		if (sample_[s].base[sorted_[2]]<sample_[s].base[sorted_[3]])
			std::swap(sorted_[2], sorted_[3]);
		if (sample_[s].base[sorted_[0]]<sample_[s].base[sorted_[1]])
			std::swap(sorted_[0], sorted_[1]);
		if (sample_[s].base[sorted_[1]]<sample_[s].base[sorted_[2]])
			std::swap(sorted_[1], sorted_[2]);
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted to access a non-existent sample_." << std::endl;
		exit(0);
	}
}

void locus::swap(const gt_t &lhs, const gt_t &rhs)
{	
	std::swap(sorted_[lhs], sorted_[rhs]);
}

count_t locus::get_count(const size_t &s, const gt_t &a) const
{
	if (a<5) return sample_[s].base[sorted_[a]];
	else return 0;
}

gt_t locus::get_index(const gt_t &a) const
{
	if (a<5) return  sorted_[a];
	else return -1;
}

count_t locus::get_count(const gt_t &a) const
{
	count_t total=0;
	if (a<5){
		for (size_t s=0; s<sample_.size() ;++s){
			if (sample_[s].masked) continue;
			total+=sample_[s].base[sorted_[a]];
		}
		return total;
	}
	else return 0;
}

count_t locus::get_coverage(const size_t &s) const
{
	if (s<sample_.size() ) return sample_[s].base[0]+
				sample_[s].base[1]+
				sample_[s].base[2]+
				sample_[s].base[3];
	else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted to access a non-existent sample_." << std::endl;
		exit(0);
	};
}

count_t locus::get_coverage() const
{
	count_t total=0;
	for (size_t s=0; s<sample_.size();++s){
		if (sample_[s].masked) continue;
		total+=sample_[s].base[0]+
			sample_[s].base[1]+
			sample_[s].base[2]+
			sample_[s].base[3];
	};
	return total;
}

/// sets quartet_t c to q.
void locus::set_quartet(const quartet_t &q, const size_t &c)
{
	sample_[c]=q;
}		

/// returns quartet_t c.
const quartet_t & locus::get_quartet(const size_t &c) const 
{
	return sample_[c];
}

/// returns quartet_t c.
quartet_t & locus::get_quartet(const size_t &c) 
{
	return sample_[c];
}
	
/// no return type.
void locus::resize(const size_t &c)
{
	sample_.resize(c);
}

/// no return type. TODO User iterator?
void locus::mask_low_cov( const count_t &dp )
{ 
	for (size_t s=0; s<sample_.size();++s) sample_[s].masked=(count(sample_[s])<dp);
}
