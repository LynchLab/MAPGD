#include "rml.h"
#include <Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

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

std::map<PAIRGL, size_t> counts;
likelihood_eq ll(0, 1);

using namespace std;

void fullModel_py (ll_t e, ll_t FA, ll_t FC, ll_t r, ll_t sA, ll_t sC, ll_t z1, ll_t z2){
	ll.fullModel(e, FA, FC, r, sA, sC, z1, z2);
}

void likelihood_eq::fullModel (ll_t p1, ll_t p2, ll_t p3, ll_t p4, ll_t p5, ll_t p6, ll_t p7, ll_t p8){
	e=p1; FA=p2; FC=p3; r=p4; sA=p5; sC=p6; z1=p7; z2=p8;
};

double get_ll_py(void){
	if (counts.size()==0) return 0;
	ll.get_ll(counts);
	return ll.ll;
};

void get_dS_py(void){
	if (counts.size()!=0) ll.get_dS(counts);
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

void likelihood_eq::get_dS(map <PAIRGL, size_t> counts){
	memset(dSll, 0, 9*sizeof(ll_t) );
	int rc;
	int i, s;
	pthread_t threads[NUM_THREADS];
	thdata data[NUM_THREADS];
	bool open[NUM_THREADS];
	pthread_attr_t attr;
	void *status;

	map<PAIRGL, size_t>::iterator start=counts.begin();
	map<PAIRGL, size_t>::iterator end=counts.end();

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
	for( i=0; i < NUM_THREADS; ++i ){
		data[i].start=start;
		data[i].eq=this;
		memset(data[i].dret, 0, 9*sizeof(ll_t) );
		open[i]=true;
		for(s=0; s < STACK_SIZE; ++s ) {++start; if(start==end) break; }
		data[i].stop=start;
		rc = pthread_create(&threads[i], NULL, &likelihood_eq::slice_dSinc, &data[i] );
		if (rc){
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}
	bool running=true;
	while(running){
		running=false;
		for( i=0; i < NUM_THREADS; i++ ){
			if (open[i]){
				rc = pthread_join(threads[i], &status);
				if (rc){
 					cout << "Error:unable to join," << rc << endl;
 					exit(-1);
 				}
 				arrayinc(dSll, data[i].dret);
 				if (start!=end){
					running=true;
					data[i].start=start;
					for(s=0; s < STACK_SIZE; ++s ) {++start; if(start==end) break; }
					data[i].stop=start;
					rc = pthread_create(&threads[i], NULL, &likelihood_eq::slice_dSinc, &data[i] );
 					pthread_attr_init(&attr);
 					pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
					if (rc){
						cout << "Error:unable to create thread," << rc << endl;
 						exit(-1);
 					}
				}
				else{
					open[i]=false;
				}
			}
		}
	}
}

void likelihood_eq::get_ll(map <PAIRGL, size_t> counts){
	ll=0;
	int rc;
	int i, s;
	pthread_t threads[NUM_THREADS];
	thdata data[NUM_THREADS];
	bool open[NUM_THREADS];
	pthread_attr_t attr;
	void *status;

	map<PAIRGL, size_t>::iterator start=counts.begin();
	map<PAIRGL, size_t>::iterator end=counts.end();

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
 
	for( i=0; i < NUM_THREADS; ++i ){
		data[i].start=start;
		data[i].eq=this;
		open[i]=true;
		for(s=0; s < STACK_SIZE; ++s ) {++start; if(start==end) break; }
		data[i].stop=start;
		rc = pthread_create(&threads[i], NULL, &likelihood_eq::slice_inc, &data[i] );
		if (rc){
			cout << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}
 
	bool running=true;
	while(running){
		running=false;
		for( i=0; i < NUM_THREADS; i++ ){
			if (open[i]){
				rc = pthread_join(threads[i], &status);
				if (rc){
 					cout << "Error:unable to join," << rc << endl;
 					exit(-1);
 				}
 				ll+=data[i].ret;
 				if (start!=end){
					running=true;
					data[i].start=start;
					for(s=0; s < STACK_SIZE; ++s ) {++start; if(start==end) break; }
					data[i].stop=start;
					rc = pthread_create(&threads[i], NULL, &likelihood_eq::slice_inc, &data[i] );
 					pthread_attr_init(&attr);
 					pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
					if (rc){
						cout << "Error:unable to create thread," << rc << endl;
 						exit(-1);
 					}
				}
				else{
					open[i]=false;
				}
			}
		}
	}
}

void* likelihood_eq::slice_inc(void *t){
	thdata* data=(thdata *) t;
	data->ret=0;
	data->it=data->start;
	while(data->it!=data->stop){
		data->ret+=data->eq->inc(data->it->first, data->it->second);
		++(data->it);
	}
	return 0;
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

	if (mm1mm2<0) mm1mm2=0;
	if (Mm1mm2<0) Mm1mm2=0;
	if (MM1mm2<0) MM1mm2=0;

	if (mm1Mm2<0) mm1Mm2=0;
	if (Mm1Mm2<0) Mm1Mm2=0;
	if (MM1Mm2<0) MM1Mm2=0;

	if (mm1MM2<0) mm1MM2=0;
	if (Mm1MM2<0) Mm1MM2=0;
	if (MM1MM2<0) MM1MM2=0;

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

int read_py (char *filename, size_t A, size_t B){
	counts.clear();
	ll.a=A;
	ll.b=B;
	counts=readcounts(filename, A, B, 0);
	return 0;
};

int read_small_py (char *filename, size_t A, size_t B, size_t s){
	counts.clear();
	ll.a=A;
	ll.b=B;
	counts=readcounts(filename, A, B, s);
	return 0;
};

ll_t dS1_py(void){return ll.dSll[0];};
ll_t dS2_py(void){return ll.dSll[1];};
ll_t dS3_py(void){return ll.dSll[2];};
ll_t dS4_py(void){return ll.dSll[3];};
ll_t dS5_py(void){return ll.dSll[4];};
ll_t dS6_py(void){return ll.dSll[5];};
ll_t dS7_py(void){return ll.dSll[6];};
ll_t dS8_py(void){return ll.dSll[7];};
ll_t dS9_py(void){return ll.dSll[8];};

int getsize_py(const char *filename){
	ifstream file;
	file.open(filename, ios::binary);
	size_t SIZE;
	file.read((char *)&SIZE,sizeof(size_t) );
	file.close();
	return SIZE;
};

PyObject * estimate_py(void){
	ll.estimate(counts);
	return Py_BuildValue("ffffffff", ll.e, ll.FA, ll.FC, ll.r, ll.sA, ll.sC, ll.z1, ll.z2);
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
		if (it->first.P!=0) N+=it->second;
		
                ++(it);
        }
	//std::cout << FA << ", " << N << std::endl;
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
	ll_t v=(P*(1-P ) );
	if(P!=0){
		FA+=( 2.*(exp(-popgl.Mm1)/4.+exp(-popgl.MM1)-P*P) -v)/v*count;
		FC+=( 2.*(exp(-popgl.Mm2)/4.+exp(-popgl.MM2)-P*P) -v)/v*count;
	};
};

inline void likelihood_eq::inc_r(const PAIRGL popgl, const ll_t count){
	ll_t P=popgl.P;
	if(P!=0){
		r+=( (exp(-popgl.Mm1-popgl.Mm2)/4.+exp(-popgl.MM1-popgl.Mm2)/2.+exp(-popgl.Mm1-popgl.MM2)/2.+exp(-popgl.MM1-popgl.MM2) )-pow(P,2) )/(P*(1-P) )*count;
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
    def("read", read_py);
    def("read_small", read_small_py);
    def("get_ll", get_ll_py);
    def("get_dS", get_dS_py);
    def("fullModel", fullModel_py);
   
    def("dS1", dS1_py);
    def("dS2", dS2_py);
    def("dS3", dS3_py);
    def("dS4", dS4_py);
    def("dS5", dS5_py);
    def("dS6", dS6_py);
    def("dS7", dS7_py);
    def("dS8", dS8_py);
    def("dS9", dS9_py);
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

std::map<PAIRGL, size_t> readcounts(const char *filename, size_t a, size_t b, size_t s){
	std::map<PAIRGL, size_t> counts;
	ifstream file;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};


	ll_t MM, Mm, mm, AN=0, BN=0, tN=0, sampled=0;
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
	
	ll_t amin=AN/tN/2.,amax=AN*2./tN, bmin=BN/2./tN, bmax=BN*2./tN;

        file.read((char *)&SIZE,sizeof(size_t) );
	if (amax>1 && bmax>1){
	        while (!file.eof()){
	                file.read((char *)&pgl.P,sizeof(ll_t) );
	                for (size_t x=0; x<SIZE; ++x){
	                        file.read((char *)&MM,sizeof(ll_t) );
	                        file.read((char *)&Mm,sizeof(ll_t) );
	                        file.read((char *)&mm,sizeof(ll_t) );
	                        file.read((char *)&N,sizeof(size_t) );
	                        pgl.add(MM, Mm, mm, N);
	                }
	                if (pgl.gl[a].N>amin && pgl.gl[a].N<amax && pgl.gl[b].N>bmin && pgl.gl[b].N<bmax){sampled++; counts[PAIRGL (pgl.gl[a].lMM, pgl.gl[a].lMm, pgl.gl[a].lmm, pgl.gl[b].lMM, pgl.gl[b].lMm, pgl.gl[b].lmm, pgl.P)]+=1;}
	                pgl.clear(); 
	   		if (sampled==s && s!=0) break;
	        };
	}
	if(s==0) cout << "Sampled: " << sampled << ", (" << a << ":" << amin << "-" << amax << ", " << b << ":" << bmin << "-" << bmax << "), ";
	file.close();
	return counts;
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
