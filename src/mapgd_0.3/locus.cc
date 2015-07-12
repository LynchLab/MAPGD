#include "locus.h"

/*! \brief site_t is initialized to id0=0, id1=0, and set to contain 0 samples.
*/	
site_t::site_t(void){
	sample.clear();
	id0=0;
	id1=0;	
};
/*! \breif site_t has a memory leak, but I don't know where it is comming from, so I've
 *	defined the crap out of it. Still leaking though.
 */
site_t::~site_t(void){
	if (sample.size()>0) sample.clear();
};

/*! \brief site_t is initialized to id0=0, id1=0, and set to contain 0 samples.
*/	
site_t::site_t(count_t size){
	sample.assign(size, quartet() );
	id0=0;
	id1=0;	
}

site_t & site_t::operator =(const site_t& arg){
        sample=arg.sample; 	
        id0=arg.id0;
	id1=arg.id1;  
        extraid=arg.extraid;   
}

void site_t::unmask(quartet_t *q){
	q->masked=false;
}

void site_t::mask(quartet_t *q){
	q->masked=true;
}

void site_t::unmask(count_t a){
	if(a<samples_) sample[a].masked=false;
}

void site_t::unmaskall(void){
	for (unsigned int s=0; s<samples_;++s){
		sample[s].masked=false;
	}
}

void site_t::maskall(void){
	for (unsigned int s=0; s<samples_;++s){
		sample[s].masked=true;
	};
};

count_t site_t::maskedcount(void) const
{
	count_t count=0;
	for(std::vector<quartet>::const_iterator it = sample.begin(); it != sample.end(); ++it) if (it->masked) count++;
	return count;
};

char site_t::getname_gt(count_t c) const
{
	if (c<5){
		if (getcount(c)==getcount(c+1) || getcount(c)==0 ) return '*';
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

char site_t::getname(count_t c) const
{
	if (c<5){
		if (getcount(c)==0 ) return '*';
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

