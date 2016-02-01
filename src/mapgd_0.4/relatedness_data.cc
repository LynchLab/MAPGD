#include "relatedness_data.h"

const std::string Relatedness::file_name=".rel";
const std::string Relatedness::table_name="RELATEDNESS";

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

void Relatedness::set_X_name(const std::string &X)
{
	X_=X;
}

void Relatedness::set_Y_name(const std::string &Y)
{
	Y_=Y;
}

std::istream& operator >> (std::istream& in, Relatedness& x) {
	std::string line;
	std::getline(in, line);
	std::stringstream stream_line(line);	
	stream_line >> x.X_ ;
	stream_line >> x.Y_ ;
	stream_line >> x.f_X_ ;
	stream_line >> x.f_X_ll ;
	stream_line >> x.f_Y_ ;
	stream_line >> x.f_Y_ll ;
	stream_line >> x.theta_XY_ ;
	stream_line >> x.theta_XY_ll ;
	stream_line >> x.gamma_XY_ ;
	stream_line >> x.gamma_XY_ll ;
	stream_line >> x.gamma_YX_ ;
	stream_line >> x.gamma_YX_ll ;
	stream_line >> x.delta_XY_ ;
	stream_line >> x.delta_XY_ll ;
	stream_line >> x.Delta_XY_ ;
	stream_line >> x.Delta_XY_ll ;
	stream_line >> x.null_ll_ ;
	stream_line >> x.max_ll_ ;
	return in;
}

std::ostream& operator<< (std::ostream& out, const Relatedness& x) {
	out << x.X_ << x.delim;
	out << x.Y_ << x.delim;
	out << x.f_X_ << x.delim;
	out << x.f_X_ll << x.delim;
	out << x.f_Y_ << x.delim;
	out << x.f_Y_ll << x.delim;
	out << x.theta_XY_ << x.delim;
	out << x.theta_XY_ll << x.delim;
	out << x.gamma_XY_ << x.delim;
	out << x.gamma_XY_ll << x.delim;
	out << x.gamma_YX_ << x.delim;
	out << x.gamma_YX_ll << x.delim;
	out << x.delta_XY_ << x.delim;
	out << x.delta_XY_ll << x.delim;
	out << x.Delta_XY_ << x.delim;
	out << x.Delta_XY_ll << x.delim;
	out << x.null_ll_ << x.delim;
	out << x.max_ll_ << x.delim;
	return out;
}

std::string Relatedness::header(void) const {
	return "@SAMPLE_X\tSAMPLE_Y\tf_X\tf_X_ll\tf_Y\tf_Y_ll\ttheta\ttheta_ll\tgamma_XY\tgamma_XY_ll\tgamma_YX\tgamma_YX_ll\tdelta\tdelta_ll\tDelta\tDelta_ll\tnull_ll\tfit\n";
}

size_t Relatedness::size(void) const {
	return sizeof(float_t)*16+(X_.size()+Y_.size())*sizeof(char);
}

void Relatedness::clear(void)
{
	X_=""; //!< the name of the first (X) sample in the compairison.
	Y_=""; //!< the name of the second (Y) sample in the compairison.

	zero();
}
void Relatedness::zero(void)
{
	sites=0;    //!< the number of sites analyzed.

	float_t e_X_[8], e_X_ll=0;
	float_t e_Y_[8], e_Y_ll=0;
	float_t f_X_=0, f_X_ll=0;
	float_t f_Y_=0, f_Y_ll=0;
	theta_XY_=0, theta_XY_ll=0;
	gamma_XY_=0, gamma_XY_ll=0;
	gamma_YX_=0, gamma_YX_ll=0;
	delta_XY_=0, delta_XY_ll=0;
	Delta_XY_=0, Delta_XY_ll=0;
	null_ll_=0, max_ll_=0;
}
void Relatedness::zero(void);  //!< zeros statistics, but doesn't set names to empty.
