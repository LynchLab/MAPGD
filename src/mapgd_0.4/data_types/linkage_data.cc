#include "linkage_data.h"

const std::string Linkage::file_name=".lds";
const std::string Linkage::table_name="LINKAGE";
const bool Linkage::binary=false;

const Registration Linkage::registered=Registration(Linkage::table_name, Linkage::create);

Linkage::Linkage ()
{
	delim='\t';
}

//TODO--IMPLEMENT
void 
Linkage::read (std::istream& in) 
{

}

//TODO--IMPLEMENT
void 
Linkage::write(std::ostream& out) const
{
	out << id1_y_ << delim;
	out << id1_y_-this->get_abs_pos() << delim;
	out << D_ << delim;
	out << Dprime() << delim;
	out << Dsq() << delim;
	out << rsq() << delim;
	out << adj_D() << delim;
	out << adj_Dprime() << delim;
	out << adj_Dsq() << delim;
	out << adj_rsq() << delim;
	out << Ni_ << delim;
	out << ll();
//	out << std::endl;
}

const std::string 
Linkage::header(void) const 
{
	return "@SCFNAME\tPOS_X\tPOS_Y\tDIST\tBEST_D\tBEST_D'\tBEST_D2\tBEST_R2\tADJ_BEST_D\tADJ_BEST_D'\tADJ_BEST_D2\tADJ_BEST_r2\tNi\tLOGLIKE\n";
}

//TODO--IMPLEMENT
const std::string 
Linkage::sql_header(void) const 
{
	return "";
}

//TODO--IMPLEMENT
const std::string 
Linkage::sql_column_names(void) const 
{
	return "";
}

//TODO--IMPLEMENT
const std::string 
Linkage::sql_values(void) const 
{
	return "";
}

//TODO--IMPLEMENT
size_t 
Linkage::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return 0;
}

const std::string 
Linkage::get_file_name(void) const
{
	return file_name;
}

const std::string 
Linkage::get_table_name(void) const
{
	return table_name;
}

float_t
Linkage::Dprime (void) const 
{
       	if (D_ >= 0) return D_/Dmax();
	return -D_/Dmin();

/* TODO: Check this. Takahiro's original function included this. I've cut it 
 * out because constructions of this form rarely have desirable statistical 
 * properties. It is possible that Dprime is an exception.
 */
/*	if (est.best_Dprime > 1.0) {
		est.best_Dprime = 1.0;
	} else if (est.best_Dprime < -1.0) {
		est.best_Dprime = -1.0;
	}*/

}

float_t 
Linkage::Dmax (void) const
{
	if ( p_*(1.0-q_) <= (1.0-p_)*q_ ) 
		return p_*(1.0-q_);
	return (1.0-p_)*q_;
}

float_t 
Linkage::Dmin (void) const
{
	if (p_*q_ <= (1.0-p_)*(1.0-q_) ) 
		return -p_*q_;
	return -(1.0-p_)*(1.0-q_);
}

float_t 
Linkage::Dsq (void) const
{
	return powl(D_,2.0);
}

float_t 
Linkage::rsq(void) const
{
	return Dsq()/(p_*(1.0-p_)*q_*(1.0-q_) );
}

float_t 
Linkage::adj_D (void) const
{
	return Ni_/(Ni_-1.0)*D_;
}

float_t 
Linkage::adj_Dsq (void) const
{
	return Ni_/(Ni_-1.0)*Dsq();
}

float_t 
Linkage::adj_Dmax (void) const
{
	return Ni_/(Ni_-1.0)*Dmax();
}

float_t 
Linkage::adj_Dmin (void) const
{
	return Ni_/(Ni_-1.0)*Dmin();
}

float_t
Linkage::adj_Dprime(void) const
{
       	if (adj_D() >= 0) return adj_D()/adj_Dmax();
	return -adj_D()/adj_Dmin();
}

float_t 
Linkage::adj_rsq (void) const
{
	return adj_Dsq()/( p_*(1.0-p_)*q_*(1.0-q_) ) - 1.0/Ni_;
}

float_t 
Linkage::ll (void) const
{
	return 2*(fit_-null_);
}
	
void Linkage::set_D(const float_t &D)
{
	D_=D;
}
void Linkage::set_fit(const float_t &fit)
{
	fit_=fit;
}
void Linkage::set_null(const float_t &null)
{
	null_=null;
}
void Linkage::set_Ni(const float_t &Ni)
{
	Ni_=Ni;
}
void Linkage::set_p(const float_t &p)
{
	p_=p;
}
void Linkage::set_q(const float_t &q)
{
	q_=q;
}

const bool Linkage::get_binary(void) const
{
	return binary;
}
