//General Likelihood functions

#include "lnmultinomial.h"

/* \breif creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .	
 */
lnmultinomial::lnmultinomial (float_t* s, const count_t& size)
{
	size_=size;
	lnp_=new float_t[size_];
	float_t *end=s+size_, *it=s, *lit=lnp_;
	while (it!=end) {*lit=log(*it); ++lit; ++it;}
}

/* \breif  . . .	
 */
lnmultinomial::lnmultinomial (const count_t &size){
	size_=size;
	lnp_=new float_t[size_];
}

/* \breif . . .	
 */
lnmultinomial::lnmultinomial (void){
	size_=4;
	lnp_=new float_t[size_];
}

/* \breif  . . .	
 */
lnmultinomial::~lnmultinomial (void){
	delete [] lnp_;
	lnfact_vector.clear();
};


/* \breif  . . .	
 */
void lnmultinomial::set(void (*fn)(allele_stat const &, float_t *), allele_stat const &s)
{
	delete [] lnp_;
	lnp_=new float_t[size_];
	fn(s, lnp_);
}

/* \breif . . .	
 */
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
}

void lnmultinomial::set(float_t a, float_t b, float_t c, float_t d)
{
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
}

float_t lnmultinomial::lnprob(const count_t *s)
{
	float_t ret=0, *lit=lnp_;
	const count_t *it=s, *end=s+size_;
	while (it!=end) {
		ret+=(*it)*(*lit); 
		++lit; 
		++it;
	}
	return lnmultinomcoef(s)+ret;
}

//cacluates the probability of a sample 
float_t lnmultinomial::lnprob_approx(const count_t *s)
{
	float_t ret=0, *lit=lnp_;
	const count_t *it=s, *end=s+size_;
	while (it!=end) {ret+=(*it)*(*lit); ++lit; ++it;}
	return ret;
}

/*string hash(count_t a){
	char *h=char *(&a);
	return *h+*(h+1);
};*/


//**NOT MULTITHREAD SAFE!!!!!!!!**/
float_t lnmultinomial::lnmultinomcoef(const count_t *s){
//**NOT MULTITHREAD SAFE!!!!!!!!**/
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

