#include "binomial.h"

uint32_t binomial::binom_coef(const uint32_t &n, const uint32_t &k){
	return fact(n)/(fact(k)*fact(n-k) );
}

uint32_t binomial::fact(const uint32_t &s){
	if (fact_vector.size()==0) { fact_vector.push_back(1); fact_vector.push_back(1);}
	for (size_t x=fact_vector.size(); x<=s; x++) fact_vector.push_back(fact_vector[x-1]*x ); 
	return fact_vector[s];
}
