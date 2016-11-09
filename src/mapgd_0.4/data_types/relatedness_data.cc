#include "relatedness_data.h"

const std::string Relatedness::file_name=".rel";
const std::string Relatedness::table_name="RELATEDNESS";
const bool Relatedness::binary=false;

const Registration Relatedness::registered=Registration(Relatedness::table_name, Relatedness::create);

Relatedness::Relatedness (){
	X_="";
	Y_="";
	delim='\t';
}

Relatedness::Relatedness (const std::string &X, const std::string &Y){
	X_=X;
	Y_=Y;
	delim='\t';
}

void 
Relatedness::set_X_name(const std::string &X)
{
	X_=X;
}

void 
Relatedness::set_Y_name(const std::string &Y)
{
	Y_=Y;
}

void
Relatedness::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::stringstream stream_line(line);	
	stream_line >> X_ ;
	stream_line >> Y_ ;
	stream_line >> f_X_ ;
	stream_line >> f_X_ll ;
	stream_line >> f_Y_ ;
	stream_line >> f_Y_ll ;
	stream_line >> theta_XY_ ;
	stream_line >> theta_XY_ll ;
	stream_line >> gamma_XY_ ;
	stream_line >> gamma_XY_ll ;
	stream_line >> gamma_YX_ ;
	stream_line >> gamma_YX_ll ;
	stream_line >> delta_XY_ ;
	stream_line >> delta_XY_ll ;
	stream_line >> Delta_XY_ ;
	stream_line >> Delta_XY_ll ;
	stream_line >> null_ll_ ;
	stream_line >> max_ll_ ;
}

void
Relatedness::write (std::ostream& out) const 
{
	out << X_ << delim;
	out << Y_ << delim;
	out << f_X_ << delim;
	out << f_X_ll << delim;
	out << f_Y_ << delim;
	out << f_Y_ll << delim;
	out << theta_XY_ << delim;
	out << theta_XY_ll << delim;
	out << gamma_XY_ << delim;
	out << gamma_XY_ll << delim;
	out << gamma_YX_ << delim;
	out << gamma_YX_ll << delim;
	out << delta_XY_ << delim;
	out << delta_XY_ll << delim;
	out << Delta_XY_ << delim;
	out << Delta_XY_ll << delim;
	out << null_ll_ << delim;
	out << max_ll_ << delim;
}

std::string 
Relatedness::header(void) const {
	return "@SAMPLE_X\tSAMPLE_Y\tf_X\tf_X_ll\tf_Y\tf_Y_ll\tθ_XY\tθ_ll\tγ_XY\tγ_XY_ll\tγ_YX\tγ_YX_ll\tδ\tδ_ll\tΔ\tΔ_ll\tnull_ll\tfit\n";
}

size_t 
Relatedness::size(void) const {
	return sizeof(float_t)*16+(X_.size()+Y_.size())*sizeof(char);
}

void 
Relatedness::clear(void)
{
	X_=""; //!< the name of the first (X) sample in the compairison.
	Y_=""; //!< the name of the second (Y) sample in the compairison.

	zero();
}

void 
Relatedness::zero(void)
{
	sites=0;    //!< the number of sites analyzed.
	memset(e_X_, 0, E_LIM*sizeof(float_t) );
	memset(e_Y_, 0, E_LIM*sizeof(float_t) );
	e_X_ll=0;
	e_Y_ll=0;
	f_X_=0, f_X_ll=0;
	f_Y_=0, f_Y_ll=0;
	theta_XY_=0, theta_XY_ll=0;
	gamma_XY_=0, gamma_XY_ll=0;
	gamma_YX_=0, gamma_YX_ll=0;
	delta_XY_=0, delta_XY_ll=0;
	Delta_XY_=0, Delta_XY_ll=0;
	null_ll_=0, max_ll_=0;
}

Relatedness& 
Relatedness::operator=(const Relatedness &rhs)
{
	sites=rhs.sites;
	memcpy(e_X_, rhs.e_X_, E_LIM*sizeof(float_t) );
	memcpy(e_Y_, rhs.e_Y_, E_LIM*sizeof(float_t) );
	e_X_ll=rhs.e_X_ll;
	e_Y_ll=rhs.e_Y_ll;
	f_X_=rhs.f_X_, f_X_ll=rhs.f_X_ll;
	f_Y_=rhs.f_Y_, f_Y_ll=rhs.f_Y_ll;
	theta_XY_=rhs.theta_XY_, theta_XY_ll=rhs.theta_XY_ll;
	gamma_XY_=rhs.gamma_XY_, gamma_XY_ll=rhs.gamma_XY_ll;
	gamma_YX_=rhs.gamma_YX_, gamma_YX_ll=rhs.gamma_YX_ll;
	delta_XY_=rhs.delta_XY_, delta_XY_ll=rhs.delta_XY_ll;
	Delta_XY_=rhs.Delta_XY_, Delta_XY_ll=rhs.Delta_XY_ll;
	null_ll_=rhs.null_ll_, max_ll_=rhs.max_ll_;
	return *this;
}

const bool 
Relatedness::get_binary(void) const
{
	return binary;
}
