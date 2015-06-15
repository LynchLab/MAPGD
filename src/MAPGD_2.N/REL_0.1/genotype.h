#ifndef _GENOTYPE_H_
#define _GENOTYPE_H_ 

#include <list>
#include <vector>
#include <tuple>
#include <iostream>

using std::size_t;

/** @brief our abstraction of log genotype probabilities.
  *
  */

class genotypelikelihood{
private :
	size_t N_;
	float_t likelihoods[3];
public :
	float_t *MM, *Mm, *mm;
	genotypelikelihood(float_t MM, float_t Mm, float_t mm, size_t N);
	genotypelikelihood();
};

class popultionlikelihood{
private:
	std::vector <genotypelikelihood>::iterator it;
public:
	std::vector <genotypelikelihood> genotypelikelihoods;
	float_t P;
	POPGL();
	popultionlikelihood(const popultionlikelihood &);
	~popultionlikelihood();
	size_t size();
	void add(genotypelikelihood gl);
	void add(float_t, float_t, float_t, size_t N);
	void clear();
};

typedef std::tuple <ll_t, ll_t, ll_t, ll_t, ll_t, ll_t, ll_t> PAIRGL; 
PAIRGL convert(POPGL&, size_t, size_t, ll_t, ll_t, ll_t, ll_t);
#endif  
