#include "compute.hpp"
#include <float.h>
#include <tuple>        
#include <algorithm>    // std::random_shuffle

using namespace std;

#define UNREL		0
#define TWINS		1
#define PAROFF		2
#define SIBS		3
#define SELFTWINS	4
#define SELFSIBS	5
#define _1COUS		6
#define _2COUS		7

likelihood_eq::likelihood_eq(size_t _a, size_t _b){
	a=_a;
	b=_b;
	ll=0;
	lastP=0;
	S1=0.0; S2=0.0; S3=0.0; S4=0.0; S5=0.0; S6=0.0; S7=0.0; S8=0.0; S9=0.0;
	b1=0.5; b2=0.5; b3=0.5; b4=0.5; b5=0.5; b6=0.5; b7=0.5; b8=0.5; b9=-4.0;
	c1=0.0; c2=0.0; c3=0.0; c4=0.0; c5=0.0; c6=0.0; c7=0.0; c8=0.0; c9=1.0;
};
void likelihood_eq::set_known(int rel){
	switch (rel){
		case UNREL:
			S9=1.0;
			S1=0;S2=0;S3=0;S4=0;S5=0;S6=0;S7=0;S8=0;
			break;
		case TWINS:
			S9=0.0;
			S1=0;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=1.0;S8=0;
			break;
		case PAROFF:
			S9=0.0;
			S1=0;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=0.0;S8=1.0;
			break;
		case SIBS:
			S9=0.25;
			S1=0;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=0.25;S8=0.5;
			break;
		case SELFSIBS:
			S9=0.0;
			S1=0.125;S2=0.125;S3=0.25;S4=0;S5=0.25;S6=0.0;S7=0.25;S8=0.0;
			break;
		case SELFTWINS:
			S9=0.0;
			S1=0.5;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=0.5;S8=0.0;
			break;
		case _1COUS:
			S9=0.75;
			S1=0.0;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=0.25;S8=0.0;
			break;
		case _2COUS:
			//S9=0.9375;
			S1=0; S2=0.15;S3=0.25;S4=0.1;S5=0.25;S6=0.2;S7=0.05;S8=0;S9=0; 
			//S1=0.0;S2=0;S3=0;S4=0;S5=0;S6=0.0;S7=0.0625;S8=0.0;
			break;
	};
};
void likelihood_eq::set(){
	S9=1.0-S1-S2-S3-S4-S5-S6-S7-S8;
};
void likelihood_eq::setmax(){
	c1=S1;c2=S2;c3=S3;c4=S4;c5=S5;c6=S6;c7=S7;c8=S8,c9=S9;
	llmax=ll;
};

/*void likelihood_eq::get95CI(){
	S1min=get;
	S1max=;
	S2min=;
	S2max=;
	S3min=;
	S3max=;
	S4min=;
	S4max=;
	S5min=;
	S5max=;
	S6min=;
	S6max=;
	S7min=;
	S7max=;
	S8min=;
	S8max=;
	S9min=;
	S9max=;
};*/

void likelihood_eq::recallmax(){
	b1=c1;b2=c2;b3=c3;b4=c4;b5=c5;b6=c6;b7=c7;b8=c8; b9=c9;
	S1=c1;S2=c2;S3=c3;S4=c4;S5=c5;S6=c6;S7=c7;S8=c8; S9=c9;
	ll=llmax;
};


void likelihood_eq::get_ll(map<PAIRGL, size_t> counts){
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
 
	//pthread_attr_destroy(&attr);

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
 
// 	pthread_exit(NULL);
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

ostream& operator<<(ostream& os, const likelihood_eq& ll){
    os << "S1: " << ll.S1;
    os << ", S2: " << ll.S2;
    os << ", S3: " << ll.S3;
    os << ", S4: " << ll.S4;
    os << ", S5: " << ll.S5;
    os << ", S6: " << ll.S6;
    os << ", S7: " << ll.S7;
    os << ", S8: " << ll.S8;
    os << ", S9: " << ll.S9;
    os << ", F " << ll.a << ":" << ll.S1+ll.S2+ll.S3+ll.S4 << ", F " << ll.b << ":" << ll.S1+ll.S2+ll.S5+ll.S6;
    os << ", R " << ll.a << "-" << ll.b << ":" << ll.S1+(ll.S3+ll.S5+ll.S7)/2.0+ll.S8/4.0;
    os << ", lnL:" << ll.ll;
    return os;

};

/*
ostream& operator<<(ostream& os, const likelihood_eq& ll){
    os << "S1: " << ll.S1 << "(" << ll.S1min "," << ll.S1max << ")"
    os << ", S2: " << ll.S2 << "(" << ll.S2min "," << ll.S2max << ")"
    os << ", S3: " << ll.S3 << "(" << ll.S3min "," << ll.S3max << ")"
    os << ", S4: " << ll.S4 << "(" << ll.S4min "," << ll.S4max << ")"
    os << ", S5: " << ll.S5 << "(" << ll.S5min "," << ll.S5max << ")"
    os << ", S6: " << ll.S6 << "(" << ll.S6min "," << ll.S6max << ")"
    os << ", S7: " << ll.S7 << "(" << ll.S7min "," << ll.S7max << ")"
    os << ", S8: " << ll.S8 << "(" << ll.S8min "," << ll.S8max << ")"
    os << ", S9: " << ll.S9 << "(" << ll.S9min "," << ll.S9max << ")"
    os << ", F" << ll.a << ":" << ll.S1+ll.S2+ll.S3+ll.S4 << ", F" << ll.b << ":" << ll.S1+ll.S2+ll.S5+ll.S6;
    os << ", R" << ll.a << "," << ll.b << ":" << ll.S1+(ll.S3+ll.S5+ll.S7)/2.0+ll.S8/4.0;
    return os;
};*/

inline ll_t likelihood_eq::inc(const PAIRGL popgl, const ll_t count){
	ll_t P, Q, Q2, Q3, Q4, P2, P3, P4, PQ;
	P=get<6>(popgl);
	Q=(1-P);
	Q2=pow((1-P), 2.0);
	Q3=pow((1-P), 3.0);
	Q4=pow((1-P),4.0);
	P2=pow(P,2.0);
	P3=pow(P,3.0);
	P4=pow(P,4.0);
	PQ=P*Q;

	ll_t MM1=get<0>(popgl), Mm1=get<1>(popgl), mm1=get<2>(popgl);
	ll_t MM2=get<3>(popgl), Mm2=get<4>(popgl), mm2=get<5>(popgl);

	return log( ( Q*S1+Q2*(S2+S3+S5+S7)+Q3*(S4+S6+S8)+Q4*S9 )*exp(-mm1-mm2)
		+( PQ*((2*S7+S8)+4*PQ*S9) )*exp(-Mm1-Mm2)
		+( P*S1+P2*(S2+S3+S5+S7)+P3*(S4+S6+S8)+P4*S9 )*exp(-MM1-MM2)

		+( PQ*S3+P*(Q2*(S4*2+S8)+2*Q3*S9) )*exp(-mm1-Mm2)
		+( PQ*S3+Q*(P2*(S4*2+S8)+2*P3*S9) )*exp(-MM1-Mm2)

		+( PQ*S5+P*(Q2*(S6*2+S8)+2*Q3*S9) )*exp(-Mm1-mm2)
		+( PQ*S5+Q*(P2*(S6*2+S8)+2*P3*S9) )*exp(-Mm1-MM2)

		+( PQ*(S2+P*S4+Q*S6)+Q2*P2*S9 )*exp(-mm1-MM2)
		+( PQ*(S2+Q*S4+P*S6)+P2*Q2*S9 )*exp(-MM1-mm2) )*count;
};

void compute (const char *filename, size_t a, size_t b, size_t subsample){
	likelihood_eq ll(a, b);

	std::map<PAIRGL, size_t> counts=readcounts(filename, a, b, subsample);

	ll_t s=0.5, maxll=-DBL_MAX;
	if (counts.size()>0){
	while (s>0.0001){
		cout << ll << endl;
		for (ll_t x1=0;x1<3;++x1){for (ll_t x2=0;x2<3;++x2){for (ll_t x3=0;x3<3;++x3){
		for (ll_t x4=0;x4<3;++x4){for (ll_t x5=0;x5<3;++x5){for (ll_t x6=0;x6<3;++x6){
		for (ll_t x7=0;x7<3;++x7){for (ll_t x8=0;x8<3;++x8){

		ll.S1=fabs(ll.b1+(x1-1.0)*s);
		ll.S2=fabs(ll.b2+(x2-1.0)*s);
		ll.S3=fabs(ll.b3+(x3-1.0)*s);
		ll.S4=fabs(ll.b4+(x4-1.0)*s);
		ll.S5=fabs(ll.b5+(x5-1.0)*s);
		ll.S6=fabs(ll.b6+(x6-1.0)*s);
		ll.S7=fabs(ll.b7+(x7-1.0)*s);
		ll.S8=fabs(ll.b8+(x8-1.0)*s);

		ll.set();
		if (ll.S9>=0){
			ll.get_ll(counts);
			if (ll.ll>maxll){
				maxll=ll.ll;
				ll.setmax();
			};
		};
		}}
		}}}
		}}};
		ll.recallmax();
		s=s*0.75;
	};
	ll.recallmax();
	cout << ll;
	ll.set_known(UNREL);
	ll.get_ll(counts);
	cout << ", Unrealted: " << ll.ll;
	ll.set_known(PAROFF);
	ll.get_ll(counts);
	cout << ", Parrent-offspring: " << ll.ll;
	ll.set_known(SIBS);
	ll.get_ll(counts);
	cout << ", Sibs: " << ll.ll;
	ll.set_known(TWINS);
	ll.get_ll(counts);
	cout << ", Twins: " << ll.ll;
	ll.set_known(TWINS);
	ll.get_ll(counts);
	cout << ", 1cous: " << ll.ll;
	ll.set_known(_1COUS);
	ll.get_ll(counts);
	cout << ", *2cous: " << ll.ll;
	ll.set_known(_2COUS);
	ll.get_ll(counts);
	cout << ", Selfed-Sibs: " << ll.ll;
	ll.set_known(SELFTWINS);
	ll.get_ll(counts);
	cout << ", Selfed-Twins: " << ll.ll << endl;
	}
	else{
		cout << "insuficent coverage\n";
	};
};
