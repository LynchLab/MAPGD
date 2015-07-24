
#include "allele_stat.h"

allele_stat::allele_stat (void){
	id0="";
	id1=0;
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

std::ostream& operator<< (std::ostream& out, const allele_stat& x) {
	if (x.coverage>0){
		out << std::fixed << std::setprecision(0);
		out << x.coverage << x.delim;
		out << std::fixed << std::setprecision(4);
		out << x.freq <<  x.delim;
		out << 1.-x.freq <<  x.delim;
		out << x.error << x.delim;
		out << x.null_error << x.delim;
		out << x.f <<  x.delim;
		out << x.MM << x.delim;
		out << x.Mm << x.delim;
		out << x.mm << x.delim;
		out << x.freq*(1-x.freq) << x.delim;
		out << (x.ll-x.monoll)*2 << x.delim;
		out << (x.ll-x.hwell)*2 << x.delim;
		out << x.gof << x.delim;
		out << x.efc << x.delim;
		out << std::fixed << std::setprecision(0);
		out << x.N << x.delim;
		out << x.excluded << x.delim;
		out << std::fixed << std::setprecision(4);
		out << x.ll;
	} else {
		out << std::fixed << std::setprecision(0);
		out << x.coverage << x.delim;
		out << '*' <<  x.delim;
		out << '*' <<  x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' <<  x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
		out << '*' << x.delim;
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

