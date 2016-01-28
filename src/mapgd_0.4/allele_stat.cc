#include "allele_stat.h"

const std::string allele_stat::file_name=".map";//!< The destination table in the Db.
const std::string allele_stat::table_name="GENOME";//!< The destination table in the Db.

allele_stat::allele_stat (void){
	id0=0;
	id1=0;
	ref=4;
	major=4;
	minor=4;
	null_error=-FLT_MAX;
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
	excluded=0;
	delim='\t';
	coverage=0;
}
	
using namespace std;

std::istream& operator >> (std::istream& in, allele_stat& x) {
	std::string line, temp;
	char c;
	in >> x.id1;
	in >> c; x.ref=Base::ctob(c);
	in >> c; x.major=Base::ctob(c);
	in >> c; x.minor=Base::ctob(c);
	in >> x.coverage;
	in >> x.freq >> temp >> x.error >> x.null_error >> x.f;
	in >> x.MM >> x.Mm >> x.mm;
	in >> temp;
	in >> x.monoll;
	in >> x.hwell;
	in >> x.gof;
	in >> x.efc;
	in >> x.N;
	in >> x.excluded;
	in >> x.ll;
	std::getline(in, line);
	x.ref=4;
	x.monoll=x.ll-x.monoll/2.;
	x.hwell=x.ll-x.hwell/2.;
	x.f-=1./(2.*x.N-1.);
	return in;
}

std::ostream& operator<< (std::ostream& out, const allele_stat& x) {
	if (x.coverage>0){
		out << x.id1 <<  x.delim;
		out << Base::btoc(x.ref) <<  x.delim;
		out << Base::btoc(x.major) <<  x.delim;
		out << Base::btoc(x.minor) <<  x.delim;
		out << std::fixed << std::setprecision(0);
		out << x.coverage << x.delim;
		out << std::fixed << std::setprecision(4);
		out << x.freq <<  x.delim;
		out << 1.-x.freq <<  x.delim;
		out << x.error << x.delim;
		out << x.null_error << x.delim;
		out << x.f+1./(2.*x.N-1.) <<  x.delim;
		out << x.MM << x.delim;
		out << x.Mm << x.delim;
		out << x.mm << x.delim;
		out << 2.*x.freq*(1.-x.freq)*(2.*x.N/(2.*x.N-1.) ) << x.delim;
		out << std::fixed << std::setprecision(2);
		out << (x.ll-x.monoll)*2 << x.delim;
		out << (x.ll-x.hwell)*2 << x.delim;
		out << x.gof << x.delim;
		out << std::setw(2);
		out << std::fixed << std::setprecision(2);
		out << x.efc << x.delim;
		out << std::fixed << std::setprecision(0);
		out << int(x.N) << x.delim;
		out << int(x.excluded) << x.delim;
		out << std::fixed << std::setprecision(4);
		out << x.ll;
	} else {
		out << x.id1 <<  x.delim;
		out << Base::btoc(x.ref) <<  x.delim;
		out << Base::btoc(x.major) <<  x.delim;
		out << Base::btoc(x.minor) <<  x.delim;
		out << std::fixed << std::setprecision(0);
		out << std::fixed << std::setprecision(0);
		out << x.coverage << x.delim;
		out << '.' <<  x.delim;
		out << '.' <<  x.delim;
		out << '.' << x.delim;
		out << '.' << x.delim;
		out << '.' <<  x.delim;
		out << '.' << x.delim;
		out << '.' << x.delim;
		out << '.' << x.delim;
		out << '.' << x.delim;
		out << std::fixed << std::setprecision(4);
		out << 0. << x.delim;
		out << 0. << x.delim;
		out << 0. << x.delim;
		out << 0. << x.delim;
		out << std::fixed << std::setprecision(0);
		out << 0 << x.delim;
		out << 0 << x.delim;
		out << std::fixed << std::setprecision(4);
		out << 0.;
	};
	return out;
};

allele_stat & allele_stat::operator=(const allele_stat & x) {
	if (this != &x) { 
		pooled=x.pooled;
		delim=x.delim;
		id0=x.id0;
		id1=x.id1;
		excluded=x.excluded;
		freq=x.freq;
		minor=x.minor;		
		major=x.major;		
		error=x.error;
		null_error=x.null_error;
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
		excluded=x.excluded;
		delim=x.delim;
	}
	return *this;
};
	
std::string allele_stat::header(void) const {
	return "@ID0    \tID1\tREF\tMAJOR\tMINOR\tCOVERAG\tMJ_FREQ\tMN_FREQ\tERROR\tNULL_ER\tF_STAT\tMM_FREQ\tMm_FREQ\tmm_FREQ\tHETERO\tPOLY_LR\tHWE_LR\tGOF\tEF_CHRM\tIND_INC\tIND_CUT\tBEST_LL\n"; 
}

size_t allele_stat::size() const{
	return sizeof(float_t)*15+sizeof(count_t)*6+sizeof(id0_t)+sizeof(char)+sizeof(bool)+sizeof(id1_t); 
}
