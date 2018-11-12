#include "allele.h"

const std::string Allele::file_name=".map";//!< The destination table in the Db.
const std::string Allele::table_name="GENOME";//!< The destination table in the Db.
const bool Allele::binary=false;

const Registration Allele::registered=Registration(Allele::table_name, Allele::create);

Allele::Allele ( std::vector<std::string> fields)
{
	abs_pos_=0;
	ref=4;
	major=4;
	minor=4;
	null_error=-FLT_MAX;
	null_error2=-FLT_MAX;
	error=-FLT_MAX;
	f=0;
	MM=1;
	Mm=0;
	mm=0;
	h=0;
	N=0;
	monoll=0;
	hwell=0;
	freq=0;
	ll=0;
	gof=0;
	efc=0;

	bias=0;
	pbias=0;

	if (fields.size()==25)
		print_bias=true;
	else
		print_bias=false;

	excluded=0;
	delim='\t';
	coverage=0;
}

Allele::Allele (void){
	abs_pos_=0;
	ref=4;
	major=4;
	minor=4;
	null_error=-FLT_MAX;
	null_error2=-FLT_MAX;
	error=-FLT_MAX;
	f=0;
	MM=1;
	Mm=0;
	mm=0;
	h=0;
	N=0;
	monoll=0;
	hwell=0;
	freq=0;
	ll=0;
	gof=0;
	efc=0;

	bias=0;
	pbias=0;
	print_bias=false;

	excluded=0;
	delim='\t';
	coverage=0;
}
	
using namespace std;

void Allele::read(std::istream& in) {
	std::string line, temp;
	std::getline(in, line);
	std::stringstream line_stream(line);
	char c;
	line_stream >> c; ref=Base::ctob(c);
	line_stream >> c; major=Base::ctob(c);
	line_stream >> c; minor=Base::ctob(c);
	line_stream >> coverage;
	//if (coverage>0){
	line_stream >> freq >> temp >> error >> null_error >> null_error2 >> f;
	line_stream >> MM >> Mm >> mm;
	line_stream >> temp;
	line_stream >> monoll;
	line_stream >> hwell;
	line_stream >> gof;
	line_stream >> efc;
	line_stream >> N;
	line_stream >> excluded;

	if (print_bias) {
		line_stream >> bias;
		line_stream >> pbias;
	}
	line_stream >> ll;

	monoll=ll-monoll/2.;
	hwell=ll-hwell/2.;
//		f-=1./(2.*N-1.);
	//}
}

/*
write_binary (std::ostream& out) const
{
	out.write( (char *) &ref, sizeof(gt_t) );
	out.write( (char *) &major, sizeof(gt_t) );
	out.write( (char *) &minor, sizeof(gt_t) );
	out.write( (char *) &coverage, sizeof() );
	out.write( () freq, );
	out.write( () error, );
	out.write( () null_error, );
	out.write( () f, );
	out.write( () MM, );
	out.write( () Mm, );
	out.write( () mm, );
	out.write( () monoll, );
	out.write( () hwell, );
	out.write( () gof, );
	out.write( () efc, );
	out.write( () N, );
	out.write( () excluded, );
	out.write( () ll, );
}*/

/*read_binary (std::istream& in)
{
	line_stream >> c; ref=Base::ctob(c);
	line_stream >> c; major=Base::ctob(c);
	line_stream >> c; minor=Base::ctob(c);
	line_stream >> coverage;
}*/

void Allele::write (std::ostream& out) const 
{
//	std::cerr << "Writing ... \n";
	if (coverage>0){
	//	std::cerr << "cov > 0 ... \n";
		out << Base::btoc(ref) <<  delim;
		out << Base::btoc(major) <<  delim;
		out << Base::btoc(minor) <<  delim;
		out << std::fixed << std::setprecision(0);
		out << coverage << delim;
		out << std::fixed << std::setprecision(4);
		out << freq <<  delim;
		if (ref==minor){
			out << freq <<  delim;
		} else if (ref==major) {
			out << 1-freq <<  delim;
		} else {
			out << 0 <<  delim;
		}
		out << error << delim;
		out << null_error << delim;
		out << null_error2 << delim;
		out << f << delim;//+1./(2.*N-1.) <<  delim;
		out << MM << delim;
		out << Mm << delim;
		out << mm << delim;
		out << 2.*freq*(1.-freq)*(2.*N/(2.*N-1.) ) << delim;
		out << std::fixed << std::setprecision(2);
		out << (ll-monoll)*2 << delim;
		out << (ll-hwell)*2 << delim;
		out << gof << delim;
		out << std::setw(2);
		out << std::fixed << std::setprecision(2);
		out << efc << delim;
		out << std::fixed << std::setprecision(0);
		out << int(N) << delim;
		out << int(excluded) << delim;
		out << std::fixed << std::setprecision(4);
	
		if (print_bias) {
			out << bias << delim;
			out << pbias << delim;
		}
		out << ll;
	} else {
	//	std::cerr << "cov == 0 ... \n";
		out << Base::btoc(ref) <<  delim;
		out << Base::btoc(ref) <<  delim;
		out << 'N' <<  delim;
		out << std::fixed << std::setprecision(0);
		out << std::fixed << std::setprecision(0);
		out << coverage << delim;
		out << '.' <<  delim;
		out << '.' <<  delim;
		out << '.' << delim;
		out << '.' << delim;
		out << '.' <<  delim;
		out << '.' << delim;
		out << '.' << delim;
		out << '.' << delim;
		out << '.' << delim;
		out << std::fixed << std::setprecision(4);
		out << 0. << delim;
		out << std::fixed << std::setprecision(2);
		out << 0. << delim;
		out << 0. << delim;
		out << 0. << delim;
		out << 0. << delim;
		out << std::fixed << std::setprecision(0);
		out << 0 << delim;
		out << 0 << delim;
		out << std::fixed << std::setprecision(4);
		if (print_bias) {
		out << 0. << delim;
		out << 0. << delim;
		}
		out << 0.;
	};
}

Allele & Allele::operator=(const Allele & x) {
	if (this != &x) { 
		pooled=x.pooled;
		delim=x.delim;
		abs_pos_=x.abs_pos_;
		ref=x.ref;
		excluded=x.excluded;
		freq=x.freq;
		minor=x.minor;		
		major=x.major;		
		error=x.error;
		null_error=x.null_error;
		null_error2=x.null_error2;
		coverage=x.coverage;	
		f=x.f;
		MM=x.MM;
		Mm=x.Mm;
		mm=x.mm;
		h=x.h;
		N=x.N;
		monoll=x.monoll;
		hwell=x.hwell;
		ll=x.ll;
		gof=x.gof;
		efc=x.efc;

		bias=x.bias;
		pbias=x.pbias;
		print_bias=x.print_bias;
	}
	return *this;
}
	
std::string Allele::header(void) const {
	if (print_bias)
		return "@SCFNAME    \tPOS\tREF\tMAJOR\tMINOR\tCOVERAG\tMJ_FREQ\tVR_FREQ\tERROR\tNULL_ER\tNULL_E2\tF_STAT\tMM_FREQ\tMm_FREQ\tmm_FREQ\tHETERO\tPOLY_LR\tHWE_LR\tGOF\tEF_CHRM\tIND_INC\tIND_CUT\tMJ_BIAS\tP_MJ_BS\tBEST_LL\n"; 
	else
		return "@SCFNAME    \tPOS\tREF\tMAJOR\tMINOR\tCOVERAG\tMJ_FREQ\tVR_FREQ\tERROR\tNULL_ER\tNULL_E2\tF_STAT\tMM_FREQ\tMm_FREQ\tmm_FREQ\tHETERO\tPOLY_LR\tHWE_LR\tGOF\tEF_CHRM\tIND_INC\tIND_CUT\tBEST_LL\n"; 
}

size_t Allele::size() const{
	return sizeof(float_t)*15+sizeof(count_t)*6+sizeof(id0_t)+sizeof(char)+sizeof(bool)+sizeof(id1_t); 
}

const std::string Allele::sql_header(void) const {
        return "(ABS_POS int, REF int, MAJOR int, MINOR int, COVERAG int, HOM_FREQ REAL, HET_FREQ REAL, POLY_LR REAL, HWE_LR REAL, GOF REAL, EF_CHRM REAL, IND_INC int, IND_CUT int, BEST_LL REAL, PRIMARY KEY (ABS_POS) )";
}

const std::string Allele::sql_column_names(void) const {
        return "(ABS_POS, REF, MAJOR, MINOR, COVERAG, HOM_FREQ, HET_FREQ, POLY_LR, HWE_LR, GOF, EF_CHRM, IND_INC, IND_CUT, BEST_LL)";
}

const std::string Allele::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
#if FLT_EVAL_METHOD == 2
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %d, %d, %d, %d, %Lf, %Lf, %Lf, %Lf, %Lf, %Lf, %d, %d, %Lf)",
	abs_pos_,
	ref,
	major,
	minor,
	coverage,
	MM,
	Mm,
	(ll-monoll)*2,
	(ll-hwell)*2,
	gof,
	efc,
	N,
	excluded,
	ll);
#elif FLT_EVAL_METHOD == 1
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %d, %d, %f)",
	abs_pos_,
	ref,
	major,
	minor,
	coverage,
	MM,
	Mm,
	(ll-monoll)*2,
	(ll-hwell)*2,
	gof,
	efc,
	N,
	excluded,
	ll);
#elif FLT_EVAL_METHOD == 0
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %d, %d, %f)",
	abs_pos_,
	ref,
	major,
	minor,
	coverage,
	MM,
	Mm,
	(ll-monoll)*2,
	(ll-hwell)*2,
	gof,
	efc,
	N,
	excluded,
	ll);
#endif
        return std::string(return_buffer);
}

const std::string Allele::get_file_name(void) const
{
	return file_name;
}

const std::string Allele::get_table_name(void) const
{
	return table_name;
}

const  bool Allele::get_binary(void) const
{
	return binary;
}
