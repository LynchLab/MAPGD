#include "genotype_pair.h"

Genotype_pair_tuple 
convert(const Genotype &x, const Genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
	return Genotype_pair_tuple ( roundf(x.lMM * T) / T, roundf(x.lMm*T)/T, roundf(x.lmm*T)/T, roundf(y.lMM*T)/T, roundf(y.lMm*T)/T, round(y.lmm*T)/T, round(m*T)/T);
}

Genotype_pair_tuple 
downvert(const Genotype &x, const Genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
	float MIN=10;
//	std::cerr << x.lMM << std::endl;
	return Genotype_pair_tuple ( roundf( std::min(x.lMM, MIN)* T) / T, roundf( std::min(x.lMm, MIN)*T)/T, roundf( std::min(x.lmm, MIN)*T)/T, roundf( std::min(y.lMM, MIN)*T)/T, roundf( std::min(y.lMm, MIN)*T)/T, round( std::min(y.lmm, MIN) *T)/T, round(m*T)/T);
}


Genotype_pair::Genotype_pair(const float_t &X_MM_, const float_t &X_Mm_, const float_t &X_mm_, const float_t &Y_MM_, const float_t &Y_Mm_, const float_t &Y_mm_, const float_t &m_)
{
	X_MM=X_MM_;
	X_Mm=X_Mm_;
	X_mm=X_mm_;
	Y_MM=Y_MM_;
	Y_Mm=Y_Mm_;
	Y_mm=Y_mm_;
	m=m_;

}

Genotype_pair_tuple 
Genotype_pair::to_tuple(const Genotype_pair &pair)
{
	return Genotype_pair_tuple (pair.X_MM, pair.X_Mm, pair.X_mm, pair.Y_MM, pair.Y_Mm, pair.Y_mm, pair.m);
}

Genotype_pair 
Genotype_pair::from_tuple(const Genotype_pair_tuple &t)
{
	return Genotype_pair (std::get<0>(t),std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), std::get<5>(t), std::get<6>(t) );
}
