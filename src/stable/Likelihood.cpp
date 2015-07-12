//General Likelihood functions

#include "Likelihood.h"
#include <iostream>

lnmultinomial::lnmultinomial (float_t* s, const count_t& size)
{
	size_=size;
	lnp_=new float_t[size_];
	float_t *end=s+size_, *it=s, *lit=lnp_;
	while (it!=end) {*lit=log(*it); ++lit; ++it;}
}

lnmultinomial::lnmultinomial (const count_t &size){
	size_=size;
	lnp_=new float_t[size_];
}

lnmultinomial::lnmultinomial (void){
	size_=4;
	lnp_=new float_t[size_];
}

lnmultinomial::~lnmultinomial (void){
	delete [] lnp_;
	lnfact_vector.clear();
};


void lnmultinomial::set(void (*fn)(allele_stat const &, float_t *), allele_stat const &s){
	delete [] lnp_;
	lnp_=new float_t[size_];
	fn(s, lnp_);
};

void lnmultinomial::set(float_t *s){
	delete [] lnp_;
	lnp_=new float_t[size_];
	float_t *end=s+size_, *it=s, *lit=lnp_;
	while (it!=end) {
		if (*it!=0) *lit=log(*it);
		else *lit=-FLT_MAX;
		++lit; 
		++it;
	}
};

void lnmultinomial::set(float_t a, float_t b, float_t c, float_t d){
	delete [] lnp_;
	lnp_=new float_t[4];
	if (a!=0) lnp_[0]=log(a);
	else lnp_[0]=-FLT_MAX;
	if (b!=0) lnp_[1]=log(b);
	else lnp_[1]=-FLT_MAX;
	if (c!=0) lnp_[2]=log(c);
	else lnp_[2]=-FLT_MAX;
	if (d!=0) lnp_[3]=log(d);
	else lnp_[3]=-FLT_MAX;
	size_=4;
};

float_t lnmultinomial::lnprob(const count_t *s){
	float_t ret=0, *lit=lnp_;
	const count_t *it=s, *end=s+size_;
	while (it!=end) {
		
		ret+=(*it)*(*lit); 
		++lit; 
		++it;
	}
	return lnmultinomcoef(s)+ret;
};

//cacluates the probability of a sample 
float_t lnmultinomial::lnprob_approx(const count_t *s){
	float_t ret=0, *lit=lnp_;
	const count_t *it=s, *end=s+size_;
	while (it!=end) {ret+=(*it)*(*lit); ++lit; ++it;}
	return ret;
};

/*string hash(count_t a){
	char *h=char *(&a);
	return *h+*(h+1);
};*/

float_t lnmultinomial::lnmultinomcoef(const count_t *s){
	return lnfact(s[0]+s[1]+s[2]+s[3])-lnfact(s[0])-lnfact(s[1])-lnfact(s[2])-lnfact(s[3]);
}

//**NOT MULTITHREAD SAFE!!!!!!!!**/
float_t lnmultinomial::lnfact(const count_t &s){
//**NOT MULTITHREAD SAFE!!!!!!!!**/
	if (lnfact_vector.size()==0) {lnfact_vector.push_back(log(1) ); lnfact_vector.push_back(log(1) );}
	for (int x=lnfact_vector.size(); x<=s; x++) lnfact_vector.push_back(lnfact_vector[x-1]+log(x) ); 
	return lnfact_vector[s];
}

struct sort_second {
		bool operator()(const std::pair<count_t, float_t> &left, const std::pair<count_t,float_t> &right) { return left.second > right.second; }
		bool operator()(const std::pair<count_t, count_t> &left, const std::pair<count_t,count_t> &right) { return left.second > right.second; }
};


std::vector <std::pair <count_t, float_t> > sort (const float_t *a, const count_t &n){
	std::vector <std::pair <count_t, float_t> > sorted;
	for (int x=0; x<n; ++x) { sorted.push_back( std::pair <count_t, float_t> (x, a[x]) ); };
	std::sort(sorted.begin(), sorted.end(), sort_second() );
	return sorted;
};

std::vector <std::pair <count_t, count_t> > sort (const count_t *a, const count_t &n){
	std::vector <std::pair <count_t, count_t> > sorted;
	for (int x=0; x<n; ++x) { sorted.push_back( std::pair <count_t, count_t> (x, a[x]) ); };
	std::sort(sorted.begin(), sorted.end(), sort_second() );
	return sorted;
}

allele_stat::allele_stat (void){
	id0="";
	id1=0;
	null_error=-FLT_MAX;
	error=-FLT_MAX;
	f=0;
	MM=1;
	Mm=0;
	mm=0;
	h=0;
	N=0;
	monoll=0;
	hwell=0;
	ll=0;
	gof=0;
	efc=0;
	excluded=0;
	delim='\t';
	coverage=0;
}
using namespace std;
std::ostream& operator<< (std::ostream& out, const allele_stat& x) {
//	out << x.id0 << x.delim;
//	out << x.id1 << x.delim;
//	out << "HI!\n";
	if (x.coverage>0){
		out << x.coverage << x.delim;
		out << x.freq <<  x.delim;
		out << 1.-x.freq <<  x.delim;
		out << x.error << x.delim;
		out << x.null_error << x.delim;
		out << x.f <<  x.delim;
		out << x.MM << x.delim;
		out << x.Mm << x.delim;
		out << x.mm << x.delim;
		out << x.h << x.delim;
		out << (x.ll-x.monoll)*2 << x.delim;
		out << (x.ll-x.hwell)*2 << x.delim;
		out << x.gof << x.delim;
		out << x.efc << x.delim;
		out << x.N << x.delim;
		out << x.excluded << x.delim;
		out << x.ll;
	} else {
		out << x.coverage << x.delim;
		out << '*' <<  x.delim;
		out << '*' <<  x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' <<  x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << 0;
	};
	return out;
};

allele_stat & allele_stat::operator=(const allele_stat & x) {
	if (this != &x) { 
		pooled=x.pooled;
		delim=x.delim;
		id0=x.id0;
		id1=x.id1;
		excluded=x.excluded;
		freq=x.freq;
		minor=x.minor;		
		major=x.major;		
		error=x.error;
		null_error=x.null_error;
		coverage=x.coverage;	
		f=x.f;
		MM=x.MM;
		Mm=x.Mm;
		mm=x.mm;
		h=x.h;
		N=x.N;
		monoll=x.monoll;
		hwell=x.hwell;
		ll=x.ll;
		gof=x.gof;
		efc=x.efc;
		excluded=x.excluded;
		delim=x.delim;
	}
	return *this;
};

