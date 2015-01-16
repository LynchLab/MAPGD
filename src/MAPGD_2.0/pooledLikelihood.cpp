#include "pooledLikelihood.h"

float_t multidist (float_t A, float_t B, float_t N, float_t a, float_t e){
	if (N==0) return 0;
	float_t pa=log( a*(1.-e)+e/3.*(1.-a) );
	float_t pb=log( (1.-a)*(1.-e)+e/3.*a );
	float_t pe=log( 2.*e/3. );
	return pa*A+pb*B+pe*(N-A-B);
};

