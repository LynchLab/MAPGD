#include "binomial.h"
uint32_t binomial::binom_coef(const uint32_t &n, const uint32_t &k){
//**NOT MULTITHREAD SAFE!!!!!!!!**/
	return fact(n)/(fact(k)*fact(n-k) );
}

float_t binomial::fact(const uint32_t &s){
	if (fact_vector.size()==0) {lnfact_vector.push_back(1); lnfact_vector.push_back(1);}
	for (int x=fact_vector.size(); x<=s; x++) lnfact_vector.push_back(fact_vector[x-1]*x ); 
	return fact_vector[s];
}
