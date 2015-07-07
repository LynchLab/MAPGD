#ifndef _COMPUTE_HPP_
#define _COMPUTE_HPP_ 1
#include <math.h>
#include <map>
#include <iostream>
#include <thread>
#include "readtable.hpp"
#include "genotype.hpp"

#define NUM_THREADS 	32
#define STACK_SIZE	100

/*
};*/
using namespace std;

class likelihood_eq{
private:
	size_t a, b;
//	ll_t P, Q, Q2, Q3, Q4, P2, P3, P4, PQ;
	inline ll_t inc(PAIRGL popgl, ll_t count);
	static void* slice_inc(void *t);
public:
	ll_t S1, S2, S3, S4, S5, S6, S7, S8, S9;
	ll_t b1, b2, b3, b4, b5, b6, b7, b8, b9;
	ll_t c1, c2, c3, c4, c5, c6, c7, c8, c9;
	ll_t lastP;
	ll_t ll, llmax;
	likelihood_eq(size_t _a, size_t _b);
	void set();
	void set_known(int);
	void get_ll(map<PAIRGL,size_t>);
	void setmax();
	void recallmax();
	friend ostream& operator<<(ostream& os, const likelihood_eq& ll);
};

typedef struct stops_thdata
{
	map<PAIRGL, size_t>::iterator start;
	map<PAIRGL, size_t>::iterator stop;
	map<PAIRGL, size_t>::iterator it;
	likelihood_eq *eq;
	ll_t ret;
} thdata;

void compute(const char* filename, size_t a, size_t b, size_t subsample);
#endif
