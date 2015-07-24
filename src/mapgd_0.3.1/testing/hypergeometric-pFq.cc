#include "hypergeometric-pFq.h"
#include <iostream>

static inline bool is_odd_D(int x) { return x % 2; }      	// note the bool

static inline int binom(int x, int n){
	while ();
	while ();
	return power_series_[x][n];
}
power_series(int x, int n){
	while ();
	while ();
	return power_series_[x][n];
};

static inline fact_sq(int n){
	while (fact_sq_.size()<=n) fact_sq_.push_back( (*fact_sq_::back)*pow(fact_sq_.size(), 2) );
	return fact_sq_[n];
};

/*! \breif for our case a1 always equals a2 and b1 is always 1, so ... */
float_t pFq::_2F1_simp(const int &a1, const float_t &z){
	std::cout << "z: " << z << std::endl;
	float_t na1=-a1;
	float_t ret=0;
	for (int n=0; n<na1; ++n){
		if is_odd_D(n) {
			ret-=(float_t)(binom(na1, n)*power_series(a1, n)*fact_sq(n) )*pow(z, n);
		} else {
			ret+=(float_t)(binom(na1, n)*power_series(a1, n)*fact_sq(n) )*pow(z, n);
		}
	}
}

float_t pFq::_2F1(const float_t &a1, const float_t &a2, const float_t &b1, const float_t &z){
	std::cout << "z: " << z << std::endl;
	return 1+a1*a2*z/b1+										//(a1_1)^2*z^2 / 1!^2
		a1*(1+a1)*a2*(1+a2)*pow(z, 2)/ (2*b1*(1+b1) )+						//(a1_2)^2*z^3 / 2!^2
		a1*(a1+1)*(a1+2)*a2*(a2+1)*(a2+2)*pow(z, 3)/(6*b1*(b1+1)*(b1+2) )+			//(a1_3)^2*z^4 / 3!^2
		a1*(a1+1)*(a1+2)*(a1+3)*a2*(a2+1)*(a2+2)*(a2+3)*pow(z, 4)/(24*b1*(b1+1)*(b1+2)*(b1+3) )	//(a1_4)^2*z^5 / 4!^2

//		a1*(a1+1)*(a1+2)*(a1+3)*a2*(a2+1)*(a2+2)*(a2+3)*pow(z, 4)/(24*b1*(b1+1)*(b1+2)*(b1+3) )	//5
//		a1*(a1+1)*(a1+2)*(a1+3)*a2*(a2+1)*(a2+2)*(a2+3)*pow(z, 4)/(24*b1*(b1+1)*(b1+2)*(b1+3) )	//5
//		a1*(a1+1)*(a1+2)*(a1+3)*a2*(a2+1)*(a2+2)*(a2+3)*pow(z, 4)/(24*b1*(b1+1)*(b1+2)*(b1+3) )	//5
		;
}

float_t pFq::_3F2(const float_t *a, const float_t *b, const float_t &z){
	return 1+a[0]*a[1]*a[2]*z/b[0]*b[1]+
		a[0]*(1+a[0])*a[1]*(1+a[1])*a[2]*(1+a[2])*pow(z, 2)/ (2*b[0]*(1+b[0])*b[1]*(1+b[1]) );
}
