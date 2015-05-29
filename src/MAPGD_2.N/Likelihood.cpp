//General Likelihood functions

#include "Likelihood.h"

#include <iostream>

lnmultinomial::lnmultinomial (float_t *s, count_t size){
	lnp=new float_t[size];
	float_t *end=s+size, *it=s, *lit=lnp;
	while (it!=end) {*lit=log(*it); ++lit; ++it;}
	_size=size;
}

lnmultinomial::lnmultinomial (count_t size){
	_size=size;
	lnp=NULL;
}

void lnmultinomial::set(void (*fn)(allele_stat const &, float_t *), allele_stat const &s){
	if (lnp!=NULL) delete lnp;
	lnp=new float_t[_size];
	fn(s, lnp);
};

void lnmultinomial::set(float_t *s){
	if (lnp!=NULL) delete lnp;
	lnp=new float_t[_size];
	float_t *end=s+_size, *it=s, *lit=lnp;
	while (it!=end) {*lit=log(*it); ++lit; ++it;}
};

void lnmultinomial::set(float_t a, float_t b, float_t c, float_t d){
	if (lnp!=NULL) delete lnp;
	lnp=new float_t[4];
	if (a!=0) lnp[0]=log(a);
	else lnp[0]=-FLT_MAX;
	if (b!=0) lnp[1]=log(b);
	else lnp[1]=-FLT_MAX;
	if (c!=0) lnp[2]=log(c);
	else lnp[2]=-FLT_MAX;
	if (d!=0) lnp[3]=log(d);
	else lnp[3]=-FLT_MAX;
	_size=4;
};

float_t lnmultinomial::lnprob(const count_t *s){
	float_t ret=0, *lit=lnp;
	const count_t *it=s, *end=s+_size;
	while (it!=end) {
		
		ret+=(*it)*(*lit); 
		++lit; 
		++it;
	}
	return lnmultinomcoef(s)+ret;
};

//cacluates the probability of a sample 
float_t lnmultinomial::lnprob_approx(const count_t *s){
	float_t ret=0, *lit=lnp;
	const count_t *it=s, *end=s+_size;
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
};


/*
allele_stat_t & allele_stat_t::operator=(const allele_stat_t& that) {
	if (this != &that) { 
		this.freq=that.freq;
		this.minor=that.minor;
		this.major=that.major;
		this.error=that.error;
		this.null_error=that.null_error;
		this.coverage=that.coverage;
		this.ll=that.ll;

		this.MM=that.MM;
		this.Mm=that.Mm;
		this.mm=that.mm;

		this.N=that.N;
		this.f=that.f;
		this.gof=that.gof;
		this.efc=that.efc;
	}
	return *this;
}*/
