#ifndef _GENOTYPE_HPP_
#define _GENOTYPE_HPP_ 

#include <list>
#include <vector>
#include <tuple>

//using std::size_t;
using namespace std;
//sitd::tuple;

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
	vector<GL>::iterator igl;
	bool frozen;
public:
	vector<GL> gl;
	ll_t P;
	POPGL();
	POPGL(const POPGL & popgl);
	~POPGL();
	size_t size();
	void add(GL gl);
	void add(ll_t lMM, ll_t lMm, ll_t lmm, size_t N);
	void clear();
};


typedef tuple <ll_t, ll_t, ll_t, ll_t, ll_t, ll_t, ll_t> PAIRGL; 
PAIRGL convert(POPGL&, size_t, size_t, ll_t, ll_t, ll_t, ll_t);
#endif  
