/* Currently likelihood calculations are all performed by a set of 'models' that we give to the "multinomial" class via the 
   multinomial::set method. These models should be able to look at the allele_stat structure, which contains information about
   the error rate at the locus and the identity of the major and minor allele, and return a set of four log probabilities of 
   observing each particular nucleotide in a given call.
*/

#include "models.h"

models& models::operator=(const models& rhs){

	if (this!=&rhs){
		lnMM_=rhs.lnMM_;
		lnMm_=rhs.lnMm_;
		lnmm_=rhs.lnmm_;

		lnMMP_=rhs.lnMMP_;
		lnMmP_=rhs.lnMmP_;
		lnmmP_=rhs.lnmmP_;
	}

	return *this;
};  

models::models(void){

	lnMM_=lnmultinomial(4);
	lnMm_=lnmultinomial(4);
	lnmm_=lnmultinomial(4);

	lnMMP_=lnmultinomial(2);
	lnMmP_=lnmultinomial(2);
	lnmmP_=lnmultinomial(2);
}

models::~models(void){
}


/*! \breif The probabilities used for calculating goodness of fit.
 *	(Major Major)
 */ 	
void MMmodelP(const allele_stat &a, float_t *l){ 	 
	if(a.error<=0){
		l[0]=0;
		l[1]=-FLT_MAX;
	} else if (a.error>=0.75) {
		l[0]=logl(0.25);
		l[1]=logl(0.75);
	} else {
		l[0]=logl(1.-a.error);
		l[1]=logl(a.error);
	}					 	
};

/*! \breif The probabilities used for calculating goodness of fit.
 *	(Major Major)
 */ 	
void mmmodelP(const allele_stat &a, float_t *l){	//Dito, assuming [m]ajor [m]inor.
	if(a.error<=0){
		l[0]=-FLT_MAX;
		l[1]=0;
	} else if (a.error>=0.75) {
		l[0]=logl(0.25);
		l[1]=logl(0.75);
	} else {
		l[0]=logl(a.error/3.);
		l[1]=logl(1-a.error/3.);
	}					 	
};

/*! \breif The probabilities used for calculating goodness of fit.
 *	(Major Major)
 */ 	
void MmmodelP(const allele_stat &a, float_t *l){ 	//[M]ajor [m]inor.
	if(a.error<=0){
		l[0]=logl(0.5);
		l[1]=logl(0.5);
	} else if (a.error>=0.75) {
		l[0]=logl(0.25);
		l[1]=logl(0.75);
	} else {
		l[0]=logl(0.5*(1.-a.error)+0.5*a.error/3.);
		l[1]=logl(0.5*(1.-a.error/3.)+0.5*a.error);
	}					 	
};

/*! \breif The probabilities used for fitting in the maximum likelihood grid search.
 *	(Major Major)
 */ 	
void MMmodel(const allele_stat &a, float_t *prob){		
	float_t e3=logl(a.error/3.);		//e3 : the error rate over three. I.e. If an error occurs, it has a 1/3 chance
						// of going to a particular base.
	if (a.error==0) e3=-FLT_MAX;		// If the error rate is zero we want to set 1/3 of the logl of error rate 
						// to the smallest (i.e. most negative) floting point number.
//	std::cout << "ln(e/3)" << e3 << std::endl;
	prob[0]=e3; prob[1]=e3; 
	prob[2]=e3; prob[3]=e3;			//set all the probabilities to the logl of the error rate over three.

	prob[a.major]=logl(1.-a.error);		//set the right error rate for major and minor.
	if (a.error==1.) prob[a.major]=-FLT_MAX;

//	std::cout << "ln(1-e)" << prob[a.major] << std::endl;
};

/*! \breif The probabilities used for fitting in the maximum likelihood grid search.
 *	(minor minor)
 */
void mmmodel(const allele_stat &a, float_t *l){	
	float_t e3=logl(a.error/3.);
	if (a.error==0.) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.minor]=logl(1.-a.error);
	if (a.error==1.) l[a.minor]=-FLT_MAX;
};

/*! \breif The probabilities used for fitting in the maximum likelihood grid search.
 *	(Major minor)
 */
void Mmmodel(const allele_stat &a, float_t *l){		//Dito.
	float_t e3=logl(a.error/3.);
	if (a.error==0.) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	float_t H=logl(0.5*(1.-a.error)+0.5*a.error/3.);
	l[a.major]=H;
	l[a.minor]=H;
	if (a.error==1.){
		l[a.major]=-FLT_MAX;
		l[a.minor]=-FLT_MAX;
	};
};

/// A function that calculates the logl likelihood of a set of observations. 
float_t models::loglikelihood(const Locus &site, const allele_stat &p){

	float_t sumll=0;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	//Tells us when to stop iterating, so we don't generate a seg. fault.

	lnMM_.set(&MMmodel, p);						//Lets initialize the multinomial distributions that will tell us the probability of  
	lnMm_.set(&Mmmodel, p);						//observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm_.set(&mmmodel, p);						//

	float_t logMM=logl(p.MM);					//The frequency of the MM genotype in the sample.
	float_t logMm=logl(p.Mm);					// Dito Mm.
	float_t logmm=logl(p.mm);					// Diot mm.

	float_t E0, E1, E2;						//These are some variables to breifly store a portion of our likelihood calculation

	while(it!=end){
		if (!it->masked){
			E0=logMM+lnMM_.lnprob_approx(it->base);
			E1=logMm+lnMm_.lnprob_approx(it->base);
			E2=logmm+lnmm_.lnprob_approx(it->base);
	
			if(E0>E2){
				if(E0>E1) sumll+=logl(expl(E1-E0)+expl(E2-E0)+1.)+E0;
				else sumll+=logl(expl(E0-E1)+expl(E2-E1)+1.)+E1;
			} 
			else if(E1>E2) sumll+=logl(expl(E0-E1)+expl(E2-E1)+1.)+E1;
			else sumll+=logl(expl(E0-E2)+expl(E1-E2)+1.)+E2;
		}
		++it;
	}
	if(isnan(sumll) ) sumll=-FLT_MAX;
	return sumll;
}

/*! \Breif DONT USE THIS!!! TODO FIX IT!. */
float_t models::genotypelikelihood(quartet_t const &quartet, const allele_stat &population){

	lnMM_.set(&MMmodel, population);	//Lets initialize the multinomial distributions that will tell us the probability of  
	lnMm_.set(&Mmmodel, population);  //observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm_.set(&mmmodel, population);  //

	float_t logMM=logl(population.MM); //The frequency of the MM genotype in the sample.
	float_t logMm=logl(population.Mm); // Dito Mm.
	float_t logmm=logl(population.mm); // Diot mm.

	float_t E[3];

	E[0]=logMM+lnMM_.lnprob(quartet.base);
	E[1]=logMm+lnMm_.lnprob(quartet.base);
	E[2]=logmm+lnmm_.lnprob(quartet.base);

	float_t ll;

	//We make E2 the largest (i.e. least negative) value of the three. E0 and E1 will often times be 
	//extreamly small, and can even be negative infinity. So we are basically just ensuring that we 
	//don't throw away a lot of precission by returning logl(exp(E2) ). 

	if(E[0]>E[2]){
		if(E[0]>E[1]) ll=logl(exp(E[1]-E[0])+expl(E[2]-E[0])+1.)+E[0];
		else ll=logl(expl(E[0]-E[1])+expl(E[2]-E[1])+1.)+E[1];
	} 
	else if(E[1]>E[2]) ll=logl(exp(E[0]-E[1])+expl(E[2]-E[1])+1.)+E[1];
	else ll=logl(expl(E[0]-E[2])+expl(E[1]-E[2])+1.)+E[2];
	
	E[0]-=ll;
	E[1]-=ll;
	E[2]-=ll;
	return E[0];
}

//THE goodness of fit calcualtion.
float_t models::goodness_of_fit (Locus &site, const allele_stat &allele, std::vector <float_t> &gofs, const float_t &MINGOF){
	float_t Num=0., Den=0., E, V, O, thisgof, clone_mingof=FLT_MAX;
	std::vector <float_t>::iterator gof=gofs.begin();

	count_t M_, N_;

	std::vector <quartet_t>::iterator it=site.sample.begin(); 
	std::vector <quartet_t>::iterator end=site.sample.end(); 

	quartet_t *maxgof_ptr=&(*it); 

	float_t logMM=logl(allele.MM); //The frequency of the MM genotype in the sample.
	float_t logMm=logl(allele.Mm); // Dito Mm.
	float_t logmm=logl(allele.mm); // Diot mm.

	lnMMP_.set(&MMmodelP, allele);  
	lnMmP_.set(&MmmodelP, allele); 
	lnmmP_.set(&mmmodelP, allele); 

	count_t count1[2];	//The count of the number of major alleles observed.

	float_t tP, etP;	//Two pointless temporary variables.

	while (it!=end){
		if (!it->masked){
			(*gof)=0;
			N_=float_t(count(*it) );
			M_=(*it).base[allele.major];

			count1[0]=M_;
			count1[1]=N_-M_;

			O=lnL(logMM, logMm, logmm, count1);

			E=0; V=0;	//The [E]xpectation and [V]arriance.

			for (size_t x=0; x<N_+1; ++x){
				count1[0]=x;			
				count1[1]=N_-x;			

				tP=lnL(logMM, logMm, logmm, count1);
				etP=exp(tP);				//

				E+=etP*tP;				//
				V+=etP*pow(tP, 2);			//
			}

			V-=pow(E, 2);
			Num+=O-E;
			Den+=V;

			if (V>0) (*gof)=( (O-E)/sqrt(V) );
			if((*gof)<clone_mingof){
				clone_mingof=(*gof);
				maxgof_ptr=&(*it);
			};
		}
		++it;
		++gof;
	}
	if (Den<=0) return 0.;
	thisgof=Num/sqrt(Den);
	if(thisgof<MINGOF){
		mask(*maxgof_ptr);
	}
	return thisgof;
}


