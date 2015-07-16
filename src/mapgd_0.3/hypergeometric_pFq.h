#ifndef HYPERGEOMETRIC_PFQ_H_
#define HYPERGEOMETRIC_PFQ_H_

#include "typedef.h"

/*! \breif pFq calculates the general hypergeometric function F. */ 
class pFq{
	private:
	public:
	static float_t _2F1(const float_t &, const float_t &, const float_t &,const float_t &);
	static float_t _3F2(const float_t *,const float_t *,const float_t &);
};
#endif
