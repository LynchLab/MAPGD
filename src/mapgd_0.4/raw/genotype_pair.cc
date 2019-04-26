#include "genotype_pair.h"

Genotype_pair_tuple 
convert(const Genotype &x, const Genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
//	std::cerr << x.MM+x.Mm*0.5 << ", " << y.MM+y.Mm*0.5 << ", " << round(m*T)/T << std::endl;
	return Genotype_pair_tuple ( roundf((x.MM) * T) / T, roundf((x.Mm)*T)/T, roundf((x.mm)*T)/T, roundf((y.MM)*T)/T, roundf((y.Mm)*T)/T, round((y.mm)*T)/T, round(m*T)/T);
}

Genotype_pair_tuple 
convert_called(const Genotype &x, const Genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
    float_t maxx=std::max({x.MM,x.Mm,x.mm});
    float_t maxy=std::max({y.MM,y.Mm,y.mm});
    return Genotype_pair_tuple ( x.MM==maxx, x.Mm==maxx, x.mm==maxx, y.MM==maxy, y.Mm==maxy, y.mm==maxy, round(m*T)/T);
}

Genotype_pair_tuple 
downvert(const Genotype &x, const Genotype &y, const float_t &m, const uint8_t &precision)
{
	float_t T=pow(10, precision);
	float_t MIN=1.2/T;
//	return Genotype_pair_tuple ( roundf( std::min((x.MM), MIN)* T) / T, roundf( std::min((x.mm), MIN)*T)/T, roundf( std::min((x.mm), MIN)*T)/T, roundf( std::min((y.MM), MIN)*T)/T, roundf( std::min((y.mm), MIN)*T)/T, round( std::min((y.mm), MIN) *T)/T, round(m*T)/T);
	return Genotype_pair_tuple ( roundf( std::max( (x.MM), MIN)* T) / T, roundf( std::max((x.Mm), MIN)*T)/T, roundf( std::max((x.mm), MIN)*T)/T, roundf( std::max((y.MM), MIN)*T)/T, roundf( std::max((y.Mm), MIN)*T)/T, round( std::max((y.mm), MIN) *T)/T, round(m*T)/T);
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

bool Genotype_pair::operator<(const Genotype_pair& rhs) const
{
	if (m!=rhs.m) return m<rhs.m;
	if (X_mm!=rhs.X_mm) return X_mm<rhs.X_mm;
	if (X_Mm!=rhs.X_Mm) return X_Mm<rhs.X_Mm;
	if (X_MM!=rhs.X_MM) return X_MM<rhs.X_MM;
	if (Y_mm!=rhs.Y_mm) return Y_mm<rhs.Y_mm;
	if (Y_Mm!=rhs.Y_Mm) return Y_Mm<rhs.Y_Mm;
	if (Y_MM!=rhs.Y_MM) return Y_MM<rhs.Y_MM;
	return false;
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
