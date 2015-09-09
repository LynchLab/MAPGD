#ifndef _RML_H_
#define _RML_H_

#include <cstring>
#include <float.h>
#include <algorithm>    // std::random_shuffle
#include <math.h>
#include <map>
#include <iostream>
#include <omp.h>
#include <list>
#include <vector>
#include <fstream>
#include <cstdio>
#include <iostream>
#include <boost/tokenizer.hpp>

#define STACK_SIZE	100

using std::size_t;

/** @brief our abstraction of log genotype probabilities.
  *
  */
typedef long double ll_t;

class GL{
private :
public :
	ll_t lMM, lMm, lmm;
	size_t N;
	GL(ll_t MM, ll_t Mm, ll_t mm, size_t lN);
	GL();
};

class POPGL{
private:
	std::vector<GL>::iterator igl;
	bool frozen;
public:
	std::vector<GL> gl;
	ll_t P;
	POPGL();
	POPGL(const POPGL & popgl);
	~POPGL();
	size_t size();
	void add(GL gl);
	void add(ll_t lMM, ll_t lMm, ll_t lmm, size_t N);
	void clear();
};

typedef struct pairgl_t{
	ll_t MM1, Mm1, mm1, MM2, Mm2, mm2, P;
	pairgl_t(ll_t MM1_, ll_t Mm1_, ll_t mm1_, ll_t MM2_, ll_t Mm2_, ll_t mm2_, ll_t P_){
		MM1=MM1_; Mm1=Mm1_; mm1=mm1_; MM2=MM2_; Mm2=Mm2_; mm2=mm2_; P=P_;
	};
	bool operator<( const pairgl_t & pg ) const {
		if (this->P==pg.P){
			if (this->MM1==pg.MM1){
				if (this->Mm1==pg.Mm1){
					if (this->mm1==pg.mm1){
						if (this->MM2==pg.MM2){
							if (this->Mm2==pg.Mm2) return this->mm2<pg.mm2;	
							else return this->Mm2<pg.Mm2;	
						}
						else return this->MM2<pg.MM2;	
					}
					else return this->mm1<pg.mm1;	
				}
				else return this->Mm1<pg.Mm1;
			}
			else return this->MM1<pg.MM1;
		}
		else return this->P<pg.P;
	}
} PAIRGL;

struct sample_stats{
	ll_t sampled, amin, amax, bmin, bmax;
	std::map <PAIRGL, size_t> counts;
};

typedef std::vector<POPGL> table_t;
table_t readtable(const char *filename);
table_t readbin(const char *filename);
sample_stats readcounts(const char *filename, size_t a, size_t b, size_t s);
void writebin(table_t table, const char *filename);

//typedef std::tuple <ll_t, ll_t, ll_t, ll_t, ll_t, ll_t> PAIRGL; 

PAIRGL convert(POPGL&, size_t, size_t, ll_t, ll_t, ll_t, ll_t);

class likelihood_eq{
private:
//	ll_t P, Q, Q2, Q3, Q4, P2, P3, P4, PQ;
	inline ll_t inc(PAIRGL popgl, ll_t count);
	inline void inc_dl(PAIRGL popgl, ll_t count);
	inline bool check(const PAIRGL &popgl);
	void inc_r(const PAIRGL popgl, const ll_t count);
	void inc_z(const PAIRGL popgl, const ll_t count);
	void inc_s(const PAIRGL popgl, const ll_t count);
	void inc_f(const PAIRGL popgl, const ll_t count);
	static void* slice_inc(void *t);
	static void* slice_dSinc(void *t);
public:
	size_t a, b, size;
	ll_t max_P;
	ll_t e, FA, FC, r, sA, sC, z1, z2;
	ll_t dfX, dfY, dr, dgXY, dgYX, dd, dD;
	ll_t J[7][7];
	ll_t H[7];

	ll_t lastP;
	ll_t ll, llmax, dSll[9];
	likelihood_eq(size_t _a, size_t _b);
	void set();
	void set_known(int);
	void get_ll(std::map<PAIRGL,size_t>);
	void get_jac(std::map<PAIRGL,size_t>);
	void get_hess(std::map<PAIRGL,size_t>);
	void fullModel(ll_t, ll_t, ll_t, ll_t, ll_t, ll_t, ll_t, ll_t);
	void setmax();
	void estimate(std::map <PAIRGL, size_t>);
	void recallmax();
	friend std::ostream& operator<<(std::ostream& os, const likelihood_eq& ll);
};

void compute(const char* filename, size_t a, size_t b, size_t subsample);

#endif
