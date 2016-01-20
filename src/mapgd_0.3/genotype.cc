#include "genotype.h"
#include <iostream>
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
