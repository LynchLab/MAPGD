#include "hypergeometric_pFq.h"

float_t pFq::_2F1(const float_t &a1, const float_t &a2, const float_t &b1, const float_t &z){
	return 1+a1*a2*z/b1+
		a1*(1+a1)*a2*(1+a2)*pow(z, 2)/ (2*b1*(1+b1) )+
		a1*(a1+1)*(a1+2)*a2*(a2+1)*(a2+2)*pow(z, 3)/(6*b1*(b1+1)*(b1+2) )
		;
}

float_t pFq::_3F2(const float_t *a, const float_t *b, const float_t &z){
	return 1+a[0]*a[1]*a[2]*z/b[0]*b[1]+
		a[0]*(1+a[0])*a[1]*(1+a[1])*a[2]*(1+a[2])*pow(z, 2)/ (2*b[0]*(1+b[0])*b[1]*(1+b[1]) );
}
