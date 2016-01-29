#include "rml.h"
#include <Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include "newton_rho.h"

using namespace std;

/** @breif default constructor. Does nothing **/
GL::GL(){};
/** @breif constuctor w/ initial values. **/
GL::GL(ll_t MM, ll_t Mm, ll_t mm, size_t lN){lMM=MM; lMm=Mm; lmm=mm; N=lN;};

POPGL::POPGL(){frozen=false;};
POPGL::POPGL(const POPGL &popgl){
	gl=popgl.gl;
	P=popgl.P;
	igl=popgl.igl;
	frozen=popgl.frozen;
};
POPGL::~POPGL(){};


likelihood_eq::likelihood_eq(size_t A, size_t B){
	a=A; b=B;	
};

/**@breif return size of POPGL if POPGL is set, 0 otherwise**/
size_t POPGL::size(){
	if (frozen) return gl.size();
	else return 0;
};

void POPGL::add(GL _gl){
	if (frozen){
		*igl=_gl;
		igl++;
	}
	else{
		gl.push_back(_gl);	
	};
};

void POPGL::add(ll_t lMM, ll_t lMm, ll_t lmm, size_t lN){
	if (frozen){
		*igl=GL(lMM, lMm, lmm, lN);
		igl++;
	}
	else{
		gl.push_back(GL(lMM, lMm, lmm, lN));	
	};
};

void POPGL::clear(){
	frozen=true;
	igl=gl.begin();
};

PAIRGL convert(POPGL& popgl, size_t a, size_t b, ll_t ma, ll_t Ma, ll_t mb, ll_t Mb){
	if ((popgl.gl[a].N>=ma) && (popgl.gl[a].N<=Ma) && (popgl.gl[b].N>=mb) && (popgl.gl[b].N<=Mb) ){
		return PAIRGL (popgl.gl[a].lMM, popgl.gl[a].lMm, popgl.gl[a].lmm, popgl.gl[b].lMM, popgl.gl[b].lMm, popgl.gl[b].lmm, popgl.P);
	}
	else return PAIRGL (0, 0, 0, 0, 0, 0, 0);
};

//Globals. counst is probably ok, since it doesn't do much, but I'm worried about likelihood_eq ll
std::map<PAIRGL, size_t> counts;
likelihood_eq ll(0, 1);

using namespace std;

void fullModel_py (ll_t e, ll_t FA, ll_t FC, ll_t r, ll_t sA, ll_t sC, ll_t z1, ll_t z2){
	ll.fullModel(e, FA, FC, r, sA, sC, z1, z2);
}

void likelihood_eq::fullModel (ll_t p1, ll_t p2, ll_t p3, ll_t p4, ll_t p5, ll_t p6, ll_t p7, ll_t p8){
	e=p1; FA=p2; FC=p3; r=p4; sA=p5; sC=p6; z1=p7; z2=p8; 
};

PyObject * get_hess_py(void){
	if (counts.size()==0) {
		std::cerr << " no data\n";
		return 0;
	}
	ll.get_hess(counts);
	return Py_BuildValue("fffffffffffffffffffffffffffffffffffffffffffffffff", 
	float(ll.J[0][0]), float(ll.J[0][1]), float(ll.J[0][2]), float(ll.J[0][3]), float(ll.J[0][4]), float(ll.J[0][5]), float(ll.J[0][6]), 
	float(ll.J[1][0]), float(ll.J[1][1]), float(ll.J[1][2]), float(ll.J[1][3]), float(ll.J[1][4]), float(ll.J[1][5]), float(ll.J[1][6]), 
	float(ll.J[2][0]), float(ll.J[2][1]), float(ll.J[2][2]), float(ll.J[2][3]), float(ll.J[2][4]), float(ll.J[2][5]), float(ll.J[2][6]), 
	float(ll.J[3][0]), float(ll.J[3][1]), float(ll.J[3][2]), float(ll.J[3][3]), float(ll.J[3][4]), float(ll.J[3][5]), float(ll.J[3][6]), 
	float(ll.J[4][0]), float(ll.J[4][1]), float(ll.J[4][2]), float(ll.J[4][3]), float(ll.J[4][4]), float(ll.J[4][5]), float(ll.J[4][6]), 
	float(ll.J[5][0]), float(ll.J[5][1]), float(ll.J[5][2]), float(ll.J[5][3]), float(ll.J[5][4]), float(ll.J[5][5]), float(ll.J[5][6]), 
	float(ll.J[6][0]), float(ll.J[6][1]), float(ll.J[6][2]), float(ll.J[6][3]), float(ll.J[6][4]), float(ll.J[6][5]), float(ll.J[6][6]) );
};

double get_max_P_py(void){
	return double(ll.max_P);
}

void set_max_P_py(ll_t max_P_py){
	ll.max_P=max_P_py;
}

PyObject * get_jac_py(void){
	if (counts.size()==0) return 0;
	ll.get_jac(counts);
	return Py_BuildValue("fffffff", float(ll.H[0]), float(ll.H[1]), float(ll.H[2]), float(ll.H[3]), float(ll.H[4]), float(ll.H[5]), float(ll.H[6]) );
};

double get_ll_py(void){
	if (counts.size()==0) return 0;
	ll.get_ll(counts);
	return double(ll.ll);
};

void arrayinc(ll_t *ret, ll_t *inc){
	ret[0]+=inc[0];
	ret[1]+=inc[1];
	ret[2]+=inc[2];
	ret[3]+=inc[3];
	ret[4]+=inc[4];
	ret[5]+=inc[5];
	ret[6]+=inc[6];
	ret[7]+=inc[7];
	ret[8]+=inc[8];
};

void likelihood_eq::get_hess(map <PAIRGL, size_t> counts){
        map<PAIRGL, size_t>::iterator start=counts.begin();
        map<PAIRGL, size_t>::iterator end=counts.end();

        vector < pair <PAIRGL, size_t> > v;

        v.reserve(counts.size() );

        for(map<PAIRGL, size_t>::iterator it=start; it!=end; ++it){
                v.push_back(*it);
        }
	
	memset(J, 0, sizeof(ll_t)*7*7);

//#	std::cout << "!!" << max_P << std::endl;

	#pragma omp parallel
	{
	ll_t J_private[7][7] = {0};
	#pragma omp for
	for (size_t z=0; z<v.size(); ++z){
		if (check(v[z].first) ) {
			for (size_t x=0; x<7; ++x){
				for (size_t y=0; y<7; ++y){
					J_private[x][y]+=Jfunc[x][y](*this, v[z].first, v[z].second);
		        	}
			}
		}
	}

	#pragma omp critical
	{
	for (size_t x=0; x<7; ++x){
		for (size_t y=0; y<7; ++y){
			J[x][y] += J_private[x][y];
		}
		}
	}
	}
}

void likelihood_eq::get_jac(map <PAIRGL, size_t> counts){

        map<PAIRGL, size_t>::iterator start=counts.begin();
        map<PAIRGL, size_t>::iterator end=counts.end();

        vector < pair <PAIRGL, size_t> > v;

        v.reserve(counts.size() );

        for(map<PAIRGL, size_t>::iterator it=start; it!=end; ++it){
                v.push_back(*it);
        }

	memset(H, 0, sizeof(ll_t)*7);
	size_t sum_out=0;

	#pragma omp parallel
	{
	ll_t H_private[7] = {0};
	size_t sum=0;
	#pragma omp for
	for (size_t x=0; x<v.size(); ++x){
		if (check(v[x].first) ) {
			for (size_t z=0; z<7; ++z){
				H_private[z]+=Hfunc[z](*this, v[x].first, v[x].second);
        		}
			sum+=v[x].second;	
		}
	}

	#pragma omp critical
	{
	for (size_t x=0; x<7; ++x){ H[x] += H_private[x]; }
	sum_out+=sum;
	}
	}
	size=sum_out;
}

void likelihood_eq::get_ll(map <PAIRGL, size_t> counts){

	map<PAIRGL, size_t>::iterator start=counts.begin();
	map<PAIRGL, size_t>::iterator end=counts.end();

	vector <pair <PAIRGL, size_t> >v;

	v.reserve(counts.size() ); 

	for(map<PAIRGL, size_t>::iterator it=start; it!=end; ++it){
		v.push_back(*it);
	} 

	ll_t sum=0;
	#pragma omp parallel for reduction(+:sum)
	for (size_t x=0; x<v.size(); ++x){
		sum+=inc(v[x].first, v[x].second);
	};
	ll=sum;
}

inline bool likelihood_eq::check(const PAIRGL &popgl){
	ll_t P, mm1mm2, Mm1mm2, MM1mm2, mm1Mm2, Mm1Mm2, MM1Mm2, mm1MM2, Mm1MM2, MM1MM2;
	P=popgl.P;

	if (max(P, 1-P) > max_P) return false;
	if (P==0 or P==1) return false;	
	
	ll_t A=P+e*P; 		//mean major allele frequency in A
	ll_t C=P-e*P; 		//mean major allele frequency in C
	ll_t Va=A*(1.-A);	//variance of the two haploid genomes of A 
	ll_t Vc=C*(1.-C);	// "   "	"	" 	"     of B
	ll_t Sa=sqrt(Va);	//standard deviation of haploid genomes A
	ll_t Sc=sqrt(Vc);	// and "	"	"	"	C.

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	ll_t E_A2  =(FA*Va+2.*pow(A,2.)+Va)/2.; //Expectation of A^2
	ll_t E_C2  =(FC*Vc+2.*pow(C,2.)+Vc)/2.; //
	ll_t E_AC  =r*Sa*Sc+A*C;
	ll_t ga=(1.-2.*A)/Sa;
	ll_t gc=(1.-2.*C)/Sc;
	ll_t E_A2C =(sA*Va*Sc*ga+A*A*C+Va*(r*(1+2.*A)+FA*C+C/(1-C) ) )/2.;
	ll_t E_AC2 =(sC*Vc*Sa*gc+C*C*A+Vc*(r*(1+2.*C)+FC*A+A/(1-A) ) )/2.;
	ll_t ka=1./(1.-A)+1./A-3.;
	ll_t kc=1./(1.-C)+1./C-3.;
	ll_t E_A2C2=(z1*sqrt(ka*kc)+z2)*Va*Vc+A*A*C*C+FA*Va*C*C+FC*Vc*A*A+4*r*Sa*Sc*A*A+C*2*sA*Va*Sc*ga+2*A*sC*Vc*Sa*gc;

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	mm1mm2=1-6*P+0.0*e+2*E_A2+2*E_C2+8.0*E_AC-4*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1mm2=0+4*P+2.0*e-4*E_A2+0*E_C2-10.*E_AC+8*E_A2C+4*E_AC2-2*E_A2C2;
	MM1mm2=0-1*P-0.5*e+2*E_A2+0*E_C2+2.0*E_AC-4*E_A2C-0*E_AC2+1*E_A2C2;
	mm1Mm2=0+4*P-2.0*e+0*E_A2-4*E_C2-10.*E_AC+4*E_A2C+8*E_AC2-2*E_A2C2;
	Mm1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2+12.*E_AC-8*E_A2C-8*E_AC2+4*E_A2C2;
	MM1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+4*E_A2C+0*E_AC2-2*E_A2C2;
	mm1MM2=0-1*P+0.5*e+0*E_A2+2*E_C2+2.0*E_AC+0*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+0*E_A2C+4*E_AC2-2*E_A2C2;
	MM1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2+0.0*E_AC+0*E_A2C+0*E_AC2+1*E_A2C2;

	if(mm1mm2>=0 and Mm1mm2>=0 and MM1mm2>=0 and mm1Mm2>=0 and Mm1Mm2>=0 and MM1Mm2>=0 and mm1MM2>=0 and Mm1MM2>=0 and MM1MM2>=0) return true;
	max_P=max(P, 1.-P);
	return false;
};

inline ll_t likelihood_eq::inc(const PAIRGL popgl, const ll_t count){
	/*This is the basic likelihood model that we are fitting. It essentially calculates the correlation coefficents
	for the first four moments of the joint distribution. The math needs to be cleaned up here. For now it is typed
	up to minimize the chance of typos. a, b, c and d are r.v. representing the presence (1) or absence of (0) of
	the Major allele in the haploid genomes of idividuals A (for a and b) and C (c and d). A is the r.v. defined as
	A~(a+b)/2 and C~(c+d)/2. Generally what we are doing here is calculating the expectations (e.g. E_A2) of the 
	joint distributions and using that to calculate the joint distribution itself. Problems can occur when 
	correlation coefficents are less than zero, becuase the probabilites of various observations (e.g. mm1mm2) can 
	become negative. These probabilities are forced to be zero, which may be a little arbitrary, but it seems to 
	work.*/


	ll_t P, mm1mm2, Mm1mm2, MM1mm2, mm1Mm2, Mm1Mm2, MM1Mm2, mm1MM2, Mm1MM2, MM1MM2;

	P=popgl.P;

	if (P==0) return 0;	
	
	ll_t A=P+e*P; 		//mean major allele frequency in A
	ll_t C=P-e*P; 		//mean major allele frequency in C
	ll_t Va=A*(1.-A);	//variance of the two haploid genomes of A 
	ll_t Vc=C*(1.-C);	// "   "	"	" 	"     of B
	ll_t Sa=sqrt(Va);	//standard deviation of haploid genomes A
	ll_t Sc=sqrt(Vc);	// and "	"	"	"	C.

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	ll_t E_A2  =(FA*Va+2.*pow(A,2.)+Va)/2.; //Expectation of A^2
	ll_t E_C2  =(FC*Vc+2.*pow(C,2.)+Vc)/2.; //
	ll_t E_AC  =r*Sa*Sc+A*C;
	ll_t ga=(1.-2.*A)/Sa;
	ll_t gc=(1.-2.*C)/Sc;
	ll_t E_A2C =(sA*Va*Sc*ga+A*A*C+Va*(r*(1+2.*A)+FA*C+C/(1-C) ) )/2.;
	ll_t E_AC2 =(sC*Vc*Sa*gc+C*C*A+Vc*(r*(1+2.*C)+FC*A+A/(1-A) ) )/2.;
	ll_t ka=1./(1.-A)+1./A-3.;
	ll_t kc=1./(1.-C)+1./C-3.;
	ll_t E_A2C2=(z1*sqrt(ka*kc)+z2)*Va*Vc+A*A*C*C+FA*Va*C*C+FC*Vc*A*A+4*r*Sa*Sc*A*A+C*2*sA*Va*Sc*ga+2*A*sC*Vc*Sa*gc;

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	mm1mm2=1-6*P+0.0*e+2*E_A2+2*E_C2+8.0*E_AC-4*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1mm2=0+4*P+2.0*e-4*E_A2+0*E_C2-10.*E_AC+8*E_A2C+4*E_AC2-2*E_A2C2;
	MM1mm2=0-1*P-0.5*e+2*E_A2+0*E_C2+2.0*E_AC-4*E_A2C-0*E_AC2+1*E_A2C2;
	mm1Mm2=0+4*P-2.0*e+0*E_A2-4*E_C2-10.*E_AC+4*E_A2C+8*E_AC2-2*E_A2C2;
	Mm1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2+12.*E_AC-8*E_A2C-8*E_AC2+4*E_A2C2;
	MM1Mm2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+4*E_A2C+0*E_AC2-2*E_A2C2;
	mm1MM2=0-1*P+0.5*e+0*E_A2+2*E_C2+2.0*E_AC+0*E_A2C-4*E_AC2+1*E_A2C2;
	Mm1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2-2.0*E_AC+0*E_A2C+4*E_AC2-2*E_A2C2;
	MM1MM2=0+0*P+0.0*e+0*E_A2+0*E_C2+0.0*E_AC+0*E_A2C+0*E_AC2+1*E_A2C2;


	if (mm1mm2<0 or mm1mm2>1) mm1mm2=0;
	if (Mm1mm2<0 or Mm1mm2>1) Mm1mm2=0; 
	if (MM1mm2<0 or MM1mm2>1) MM1mm2=0; 
	if (mm1Mm2<0 or mm1Mm2>1) mm1Mm2=0; 
	if (Mm1Mm2<0 or Mm1Mm2>1) Mm1Mm2=0;
	if (MM1Mm2<0 or MM1Mm2>1) MM1Mm2=0;
	if (mm1MM2<0 or mm1MM2>1) mm1MM2=0;
	if (Mm1MM2<0 or Mm1MM2>1) Mm1MM2=0;
	if (MM1MM2<0 or MM1MM2>1) MM1MM2=0;

	ll_t S=pow(mm1mm2+Mm1mm2+MM1mm2+mm1Mm2+Mm1Mm2+MM1Mm2+mm1MM2+Mm1MM2+MM1MM2, 2);

	if(S>1){
		mm1mm2/=S;
		Mm1mm2/=S;
		MM1mm2/=S;
		mm1Mm2/=S;
		Mm1Mm2/=S;
		MM1Mm2/=S;
		mm1MM2/=S;
		Mm1MM2/=S;
		MM1MM2/=S;
	};
	
	ll_t E[9];

        if (mm1mm2>0) E[0]=log(mm1mm2)-popgl.mm1-popgl.mm2;
     	else E[0]=-FLT_MAX;
	if (Mm1Mm2>0) E[1]=log(Mm1Mm2)-popgl.Mm1-popgl.Mm2;
     	else E[1]=-FLT_MAX;
        if (MM1MM2>0) E[2]=log(MM1MM2)-popgl.MM1-popgl.MM2;
     	else E[2]=-FLT_MAX;
        if (mm1Mm2>0) E[3]=log(mm1Mm2)-popgl.mm1-popgl.Mm2;
     	else E[3]=-FLT_MAX;
        if (MM1Mm2>0) E[4]=log(MM1Mm2)-popgl.MM1-popgl.Mm2;
     	else E[4]=-FLT_MAX;
        if (Mm1mm2>0) E[5]=log(Mm1mm2)-popgl.Mm1-popgl.mm2;
     	else E[5]=-FLT_MAX;
        if (Mm1MM2>0) E[6]=log(Mm1MM2)-popgl.Mm1-popgl.MM2;
     	else E[6]=-FLT_MAX;
        if (mm1MM2>0) E[7]=log(mm1MM2)-popgl.mm1-popgl.MM2;
     	else E[7]=-FLT_MAX;
        if (MM1mm2>0) E[8]=log(MM1mm2)-popgl.MM1-popgl.mm2;
     	else E[8]=-FLT_MAX;

	std::sort(E, E+9);

	return (log(1+exp(E[0]-E[8])+exp(E[1]-E[8])+exp(E[2]-E[8])+exp(E[3]-E[8])+exp(E[4]-E[8])+exp(E[5]-E[8])+exp(E[6]-E[8])+exp(E[7]-E[8]) )+E[8] )*count;
};

PyObject * read_py (char *filename, size_t A, size_t B){
	counts.clear();
	ll.a=A;
	ll.b=B;
	sample_stats sample=readcounts(filename, A, B, 0);
	counts=sample.counts;
	return Py_BuildValue("fffffff", float(sample.sampled), float(A), float(sample.amin), float(sample.amax), float(B), float(sample.bmin), float(sample.bmax) );
};

PyObject * read_small_py (char *filename, size_t A, size_t B, size_t s){
	counts.clear();
	ll.a=A;
	ll.b=B;
	sample_stats sample=readcounts(filename, A, B, 0);
	counts=sample.counts;
	return Py_BuildValue("fffff", float(sample.sampled), float(sample.amin), float(sample.amax), float(sample.bmin), float(sample.bmax) );
};

/*
ll_t dS1_py(void){return ll.dSll[0];};
ll_t dS2_py(void){return ll.dSll[1];};
ll_t dS3_py(void){return ll.dSll[2];};
ll_t dS4_py(void){return ll.dSll[3];};
ll_t dS5_py(void){return ll.dSll[4];};
ll_t dS6_py(void){return ll.dSll[5];};
ll_t dS7_py(void){return ll.dSll[6];};
ll_t dS8_py(void){return ll.dSll[7];};
ll_t dS9_py(void){return ll.dSll[8];};
*/

int size_py(void){
	return int(ll.size);
}

int getsize_py(const char *filename){
	ifstream file;
	file.open(filename, ios::binary);
	size_t SIZE;
	file.read((char *)&SIZE,sizeof(size_t) );
	file.close();
	return int(SIZE);
};

PyObject * estimate_py(void){
	ll.estimate(counts);
	return Py_BuildValue("ffffffff", float(ll.e), float(ll.FA), float(ll.FC), float(ll.r), float(ll.sA), float(ll.sC), float(ll.z1), float(ll.z2) );
};

void likelihood_eq::estimate(map <PAIRGL, size_t> counts){
	e=0; FA=0; FC=0; r=0; sA=0; sC=0; z1=0; z2=0;
	map<PAIRGL, size_t>::iterator start=counts.begin();
	map<PAIRGL, size_t>::iterator end=counts.end();
	map<PAIRGL, size_t>::iterator it=start;

	ll_t N=0;
        while(it!=end){
                inc_f(it->first, it->second);
                inc_r(it->first, it->second);
		if (it->first.P!=0) N+=it->second/(it->first.P*(1-it->first.P) );
		
                ++(it);
        }
	FA/=N;
	FC/=N;
	r/=N;
	N=0;
	it=start;
        while(it!=end){
                inc_s(it->first, it->second);
		if (it->first.P!=0 && it->first.P!=0.5) N+=it->second;
                ++(it);
        }
	sA/=N;
	sC/=N;
	it=start;
        while(it!=end){
                inc_z(it->first, it->second);
                ++(it);
        }
};

inline void likelihood_eq::inc_f(const PAIRGL popgl, const ll_t count){
	ll_t P=popgl.P;
	ll_t v=(P*(1.-P ) );
	if(P!=0){
		FA+=( 2.*(exp(-popgl.Mm1)/4.+exp(-popgl.MM1)-P*P) -v)/pow(v,2)*count;
		FC+=( 2.*(exp(-popgl.Mm2)/4.+exp(-popgl.MM2)-P*P) -v)/pow(v,2)*count;
	};
};

inline void likelihood_eq::inc_r(const PAIRGL popgl, const ll_t count){
	ll_t P=popgl.P;
	if(P!=0){
		r+=( (exp(-popgl.Mm1-popgl.Mm2)/4.+exp(-popgl.MM1-popgl.Mm2)/2.+exp(-popgl.Mm1-popgl.MM2)/2.+exp(-popgl.MM1-popgl.MM2) )-pow(P,2) )/pow(P*(1-P),2)*count;
	};
};

inline void likelihood_eq::inc_s(const PAIRGL popgl, const ll_t count){
	ll_t P=popgl.P;
	if(P!=0 && P!=0.5){
		sA+=( 2.*(exp(-popgl.Mm1-popgl.Mm2)/8.+exp(-popgl.MM1-popgl.Mm2)/2.+exp(-popgl.Mm1-popgl.MM2)/4.+exp(-popgl.MM1-popgl.MM2) )-pow(P,3)-P*(1-P)*(r*(1+2*P)+FA*P+P/(1-P) ) )/(P*(1-P) )/(0.5-P)/2*count;
		sC+=( 2.*(exp(-popgl.Mm1-popgl.Mm2)/8.+exp(-popgl.MM1-popgl.Mm2)/4.+exp(-popgl.Mm1-popgl.MM2)/2.+exp(-popgl.MM1-popgl.MM2) )-pow(P,3)-P*(1.-P)*(r*(1.+2.*P)+FC*P+P/(1.-P) ) )/(P*(1.-P) )/(0.5-P)/2.*count;
	};
}

inline void likelihood_eq::inc_z(const PAIRGL popgl, const ll_t count){
	ll_t P=popgl.P;
}

using namespace boost::python;

BOOST_PYTHON_MODULE(rml)
{
    def("estimate", estimate_py);
    def("getsize", getsize_py);
    def("size", size_py);
    def("read", read_py);
    def("get_ll", get_ll_py);
    def("get_max_P", get_max_P_py);
    def("set_max_P", set_max_P_py);
    def("get_jac", get_jac_py);
    def("get_hess", get_hess_py);
    def("fullModel", fullModel_py);
};

const unsigned char MIN_QUAL = '!';

ll_t max(ll_t a, ll_t b){
	if (a>b) return a;
	else return b;
}; 

table_t readtable(const char *filename){
	ifstream file;
	string buffer;
	string token;
	file.open(filename);
	table_t table;
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	//our first order of buisness is to figure out the number of coloumns in the file;
	ll_t MM, Mm, mm, P;
	size_t N;
	POPGL pgl;
	while (!file.eof()){
		getline(file,buffer);
		if (file.eof() ) break;
   		boost::tokenizer<> tok(buffer);
		boost::char_separator<char> sep("\t");
		boost::tokenizer<boost::char_separator<char> > tokens(buffer, sep);
		boost::tokenizer<boost::char_separator<char> >::iterator t=tokens.begin();
		for (int x=0; x<4; ++x) ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
		P=atof(t->c_str());
		pgl.P=P;
		for (int x=0; x<2; ++x) ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
		while (t!=tokens.end() ){
			MM=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			Mm=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			mm=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			N=atoi(t->c_str()); ++t;
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(POPGL(pgl));
		pgl.clear(); 
		buffer="";
	};
	return table;
};

table_t readbin(const char *filename){
	ifstream file;
	string buffer;
	file.open(filename, ios::binary);
	table_t table;
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	ll_t MM, Mm, mm;
	size_t N, SIZE;
	file.read((char *)&SIZE,sizeof(size_t) );
	POPGL pgl;
	while (!file.eof()){
		file.read((char *)&pgl.P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(ll_t) );
			file.read((char *)&Mm,sizeof(ll_t) );
			file.read((char *)&mm,sizeof(ll_t) );
			file.read((char *)&N,sizeof(size_t) );
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(POPGL(pgl));
		pgl.clear(); 
	};
	file.close();
	return table;
};


sample_stats readcounts (const char *filename, size_t a, size_t b, size_t s){
	sample_stats sample;
	std::map <PAIRGL, size_t> counts;
	ifstream file;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};


	ll_t MM, Mm, mm, AN=0, BN=0, tN=0;
	size_t SIZE, N;
	file.read((char *)&SIZE,sizeof(size_t) );
	POPGL pgl;
	while (!file.eof()){
		file.read((char *)&pgl.P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(ll_t) );
			file.read((char *)&Mm,sizeof(ll_t) );
			file.read((char *)&mm,sizeof(ll_t) );
			file.read((char *)&N,sizeof(size_t) );
			if (x==a) AN+=N;
			if (x==b) BN+=N;
		}
		tN+=1;
	};
	file.close();
	file.open(filename, ios::binary);
	
	sample.amin=AN/2./tN;
	sample.amax=AN*2./tN;
	sample.bmin=BN/2./tN;
	sample.bmax=BN*2./tN;

        file.read((char *)&SIZE,sizeof(size_t) );
	sample.sampled=0;

	if (sample.amax>1 && sample.bmax>1){
	        while (!file.eof()){
	                file.read((char *)&pgl.P,sizeof(ll_t) );
	                for (size_t x=0; x<SIZE; ++x){
	                        file.read((char *)&MM,sizeof(ll_t) );
	                        file.read((char *)&Mm,sizeof(ll_t) );
	                        file.read((char *)&mm,sizeof(ll_t) );
	                        file.read((char *)&N,sizeof(size_t) );
	                        pgl.add(MM, Mm, mm, N);
	                }
	                if (pgl.gl[a].N>sample.amin && pgl.gl[a].N<sample.amax && pgl.gl[b].N>sample.bmin && pgl.gl[b].N<sample.bmax){
				sample.sampled++; 
				sample.counts[PAIRGL (pgl.gl[a].lMM, pgl.gl[a].lMm, pgl.gl[a].lmm, pgl.gl[b].lMM, pgl.gl[b].lMm, pgl.gl[b].lmm, pgl.P)]+=1;
			}
	                pgl.clear(); 
	   		if (sample.sampled==s && s!=0) break;
	        };
	}
	file.close();
	return sample;
};

void writebin(table_t table, const char *filename){
	ofstream file;
	string buffer;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	size_t SIZE;
	SIZE=table[1].size();

	table_t::iterator it=table.begin();
	table_t::iterator end=table.end();
	file.write((char *)&SIZE,sizeof(size_t) );
	POPGL pgl;
	while (it!=end){
		file.write((char *)&it->P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.write((char *)&it->gl[x].lMM,sizeof(ll_t) );
			file.write((char *)&it->gl[x].lMm,sizeof(ll_t) );
			file.write((char *)&it->gl[x].lmm,sizeof(ll_t) );
			file.write((char *)&it->gl[x].N,sizeof(size_t) );
		}
		it++;
	};
};
