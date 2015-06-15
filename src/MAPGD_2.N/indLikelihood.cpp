#include "indLikelihood.h"
/* Written by Matthew Ackerman.*/

/* Currently likelihood calculations are all performed by a set of 'models' that we give to the "multinomial" class via the 
   multinomial::set method. These models should be able to look at the allele_stat structure, which contains information about
   the error rate at the locus and the identity of the major and minor allele, and return a set of four log probabilities of 
   observing each particular nucleotide in a given call.
*/

void MMmodelP(allele_stat const &a, float_t *l){ 	//The model used for the goodness of fit test assuming the genotype is 
						 	//[M]ajor [M]ajor.
	float_t H=(1-a.error);
	float_t e3=log(1-H);
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
};

void mmmodelP(allele_stat const &a, float_t *l){	//Dito, assuming [m]ajor [m]inor.
	float_t H=(a.error/3);
	float_t e3=log(1-H);
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
};

void MmmodelP(allele_stat const &a, float_t *l){ 	//[M]ajor [m]inor.
	float_t H=(0.5*(1.-a.error)+0.5*a.error/3.);
	float_t e3=log(1-H);
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
};

void MMmodel(allele_stat const &a, float_t *l){		//These are the full probabilities used for fitting the maximum
	float_t e3=log(a.error/3.);			//likelihood model.
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(1.-a.error);
};

void mmmodel(allele_stat const &a, float_t *l){		//Dito.
	float_t e3=log(a.error/3.);
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.minor]=log(1.-a.error);
};

void Mmmodel(allele_stat const &a, float_t *l){		//Dito.
	float_t e3=log(a.error/3.);
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	float_t H=log(0.5*(1.-a.error)+0.5*a.error/3.);
	l[a.major]=H;
	l[a.minor]=H;
};

lnmultinomial lnMM(4), lnMm(4), lnmm(4), lnF(4);	//The three multinomials we will use for probability calculations.
							//The '4' specifies the number of categories of the distribution.
							//Since these represent the distribution of the four nucleotides 
							//A, C, G and T, we use 4 categories.

/*@Breif, a function for calculating the jacobian of the log likelihood something something. */
float_t dll(profile const &pro, allele_stat const &p, count_t const &MIN){
}

/*@Breif, a function that calculats the log likelihood of a set of observations. */
float_t loglikelihood(profile const &pro, allele_stat const &p, count_t const &MIN){

	float_t sumll=0;

	std::vector <quartet_t>::const_iterator it=pro.begin();	//Lets us iterate of the quartets. 
	std::vector <quartet_t>::const_iterator end=pro.end();	//Tells us when to stop iterating, so we don't generate a seg. fault.

	lnMM.set(&MMmodel, p);	//Lets initialize the multinomial distributions that will tell us the probability of  
	lnMm.set(&Mmmodel, p);  //observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm.set(&mmmodel, p);  //

	float_t logMM=log(p.MM); //The frequency of the MM genotype in the sample.
	float_t logMm=log(p.Mm); // Dito Mm.
	float_t logmm=log(p.mm); // Diot mm.

	float_t E0, E1, E2; //These are some variables to breifly store a portion of our likelihood calculation

	count_t T;

	while(it!=end){
		if (!it->masked){
			T=( it->base[0]+it->base[1]+it->base[2]+it->base[3] );
			if ( T>=MIN ) {

				E0=logMM+lnMM.lnprob(it->base);
				E1=logMm+lnMm.lnprob(it->base);
				E2=logmm+lnmm.lnprob(it->base);

				//We make E2 the largest (i.e. least negative) value of the three. E0 and E1 will often times be 
				//extreamly small, and can even be negative infinity. So we are basically just ensuring that we 
				//don't throw away a lot of precission by returning log(exp(E2) ). 

				if(E0>E2){
					if(E0>E1) sumll+=log(exp(E1-E0)+exp(E2-E0)+1.)+E0;
					else sumll+=log(exp(E0-E1)+exp(E2-E1)+1.)+E1;
				} 
				else if(E1>E2) sumll+=log(exp(E0-E1)+exp(E2-E1)+1.)+E1;
				else sumll+=log(exp(E0-E2)+exp(E1-E2)+1.)+E2;
			}
		}
		++it;
	}
	return sumll;
}

/*@Breif, a function that calculats the log likelihood of a set of observations. */
float_t genotypelikelihood(quartet_t const &quartet, allele_stat const &population, count_t const &MIN){

	lnMM.set(&MMmodel, population);	//Lets initialize the multinomial distributions that will tell us the probability of  
	lnMm.set(&Mmmodel, population);  //observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm.set(&mmmodel, population);  //

	float_t logMM=log(population.MM); //The frequency of the MM genotype in the sample.
	float_t logMm=log(population.Mm); // Dito Mm.
	float_t logmm=log(population.mm); // Diot mm.

	float_t E[3];

	E[0]=logMM+lnMM.lnprob(quartet.base);
	E[1]=logMm+lnMm.lnprob(quartet.base);
	E[2]=logmm+lnmm.lnprob(quartet.base);

	float_t ll;

	//We make E2 the largest (i.e. least negative) value of the three. E0 and E1 will often times be 
	//extreamly small, and can even be negative infinity. So we are basically just ensuring that we 
	//don't throw away a lot of precission by returning log(exp(E2) ). 

	if(E[0]>E[2]){
		if(E[0]>E[1]) ll=log(exp(E[1]-E[0])+exp(E[2]-E[0])+1.)+E[0];
		else ll=log(exp(E[0]-E[1])+exp(E[2]-E[1])+1.)+E[1];
	} 
	else if(E[1]>E[2]) ll=log(exp(E[0]-E[1])+exp(E[2]-E[1])+1.)+E[1];
	else ll=log(exp(E[0]-E[2])+exp(E[1]-E[2])+1.)+E[2];
	
	E[0]-=ll;
	E[1]-=ll;
	E[2]-=ll;
}

/*@Breif, Used for GOF, same basic idea as above. I know, if I have to write it twice I'm doing something wrong.*/
float_t lnP(count_t const *base, allele_stat const &p){

	lnMM.set(&MMmodelP, p);
	lnMm.set(&MmmodelP, p);
	lnmm.set(&mmmodelP, p);
	
	float_t logMM=log(p.MM);
	float_t logMm=log(p.Mm);
	float_t logmm=log(p.mm);

	float_t E0, E1, E2, tE;

	E0=logMM+lnMM.lnprob(base);
	E1=logMm+lnMm.lnprob(base);
	E2=logmm+lnmm.lnprob(base);

	if(E0>E2) {tE=E2; E2=E0; E0=tE;}
	if(E1>E2) {tE=E2; E2=E1; E1=tE;}

	return log(exp(E0-E2)+exp(E1-E2)+1.)+E2;
};

/*Used for GOF, this calculate the expectation ad variance. The majority of the range of x may contribute essentially nothing 
  to this calculation, so I may look into some better bounds to speed up these calculations.*/

void EVofP (count_t N_, allele_stat const &a, float_t &E, float_t &V){
	count_t l[4];
	float_t tP;
	E=0; V=0;
	for (int x=0; x<N_+1; ++x){
		memset(l, 0, sizeof(count_t)*4);
		l[a.major]=x;
		l[a.minor]=N_-x;
		tP=lnP(l, a);
		E+=exp(tP)*tP;
		V+=exp(tP)*pow(tP, 2);
	}
	V-=pow(E, 2);
}

//THE goodness of fit calcualtion.
float_t gof (profile &pro, allele_stat const &a, count_t const &MIN, float_t const &site_maxgof){
	float_t Num=0., Den=0., E, V, O, thisgof, gof, clone_mingof=FLT_MAX;

	count_t M_, N_, l[4];

	std::vector <quartet_t>::iterator it=pro.begin(); 
	std::vector <quartet_t>::iterator end=pro.end(); 
	quartet_t * maxgof_ptr; 

	while (it!=end){
		if (!it->masked){
			N_=float_t(count(*it));
			if ( N_>MIN ){
				memset(l, 0, sizeof(count_t)*4);
				M_=(*it).base[a.major];
				l[a.major]=M_;
				l[a.minor]=N_-M_;
				EVofP(N_, a, E, V);
				O=lnP(l, a);
				//std::cout << M_ << ", " << N_ << ", " << O << ", " << E << ", " << V <<std::endl;
				Num+=O-E;
				Den+=V;
				gof=(O-E)/sqrt(V);
				if(gof<clone_mingof){
					clone_mingof=gof;
					maxgof_ptr=&(*it);
				};
			}
		}
		++it;
	}
	if (Den<=0) return 0.;
	thisgof=Num/sqrt(Den);
	if(abs(thisgof)>site_maxgof) pro.mask(maxgof_ptr);
	return thisgof;
}

/*calculates the 'effective number of chromosomes' in a sample. This number isn't used for anything right now.*/

float_t efc (profile const &pro, count_t const &MIN){
	std::vector <quartet_t>::const_iterator it=pro.begin(); 
	std::vector <quartet_t>::const_iterator end=pro.end(); 

	float_t ec=0;
	count_t N_;

	while (it!=end){
		if (!it->masked){
			N_=count(*it);
			if ( N_>MIN ){
				ec+=2-2*pow(0.5, N_);
			}
		}
		++it;
	}
	return ec;
};

count_t initparams_old(profile &pro, allele_stat &a, const count_t &MIN, const float_t &minerr, const count_t &M){

	std::vector <quartet_t>::const_iterator it=pro.begin(); 
	std::vector <quartet_t>::const_iterator end=pro.end(); 

	it=pro.begin(); 

	count_t allele_freq[4], genotype[4][4], N_=0;
	
	memset(allele_freq,0, 4*sizeof(count_t) );

	memset(genotype,0, 4*sizeof(count_t) );
	memset(genotype+1,0, 4*sizeof(count_t) );
	memset(genotype+2,0, 4*sizeof(count_t) );
	memset(genotype+3,0, 4*sizeof(count_t) );

	while (it!=end){
		if (!it->masked){
			float_t T=float_t(count(*it));
			if (T>MIN){
				N_+=1.;
				for (count_t x=0; x<4; ++x){
					if (float_t( (*it).base[x])/T>0.9) {allele_freq[x]+=2.; genotype[x][x]+=1.;}
					for (count_t y=x+1; y<4; ++y){
						if (float_t( (*it).base[y])/T<0.9)
						if (float_t( (*it).base[x])/T<0.9)
						if (float_t( (*it).base[y]+(*it).base[x])/T>0.9 ) {allele_freq[x]+=1.; allele_freq[y]+=1.; genotype[x][y]+=1.;}
					};
				};
			}
		}
		++it;
	}

	if (N_==0) return 0;	
	
	std::vector <std::pair <count_t, count_t> > sorted=sort(allele_freq, 4);	

	a.major=sorted[0].first;
	a.minor=sorted[M+1].first;

	count_t S=pro.getcoverage();
	count_t M_=pro.getcount(a.major);
	count_t m_=pro.getcount(a.minor);

	count_t E_=S-M_-m_;

	float_t e_;

	if (minerr>0){
		e_= float_t(E_)/float_t(S);
		a.error=3.*e_/2.;
		e_=float_t(E_+m_)/float_t(S);
		a.null_error=e_;
		if (a.error<minerr) a.error=minerr;
		if (a.null_error<minerr) a.null_error=minerr;
	} else {
		e_= float_t(E_)/float_t(S);
		a.error=3.*e_/2.;
		e_=float_t(E_+m_)/float_t(S);
		a.null_error=e_;
	};
	float_t P_=genotype[a.major][a.major], Q_=genotype[a.minor][a.minor], H_=genotype[a.major][a.minor]+genotype[a.minor][a.major];
	a.N=N_;
	N_=P_+Q_+H_;
	P_/=N_;
	Q_/=N_;
	H_/=N_;
	a.freq=P_+H_/2;
	a.MM=P_;
	a.Mm=H_;
	a.mm=Q_;
	a.efc=efc(pro, MIN);
	if (M==2) return 0;
	if (allele_freq[M+2]!=0 and allele_freq[M+1]==allele_freq[M+2]) return 1;
	return 0;
}

//Intilaizes the parameters of the ... returns 0 on succesful excecution, returns 1 if there is am...
count_t initparams(profile &pro, allele_stat &a, const count_t &MIN, const float_t &minerr, const count_t &M){

	std::vector <quartet_t>::const_iterator it=pro.begin(); 
	std::vector <quartet_t>::const_iterator end=pro.end(); 

	it=pro.begin(); 
	pro.sort();

	count_t N_=0;


	count_t S=pro.getcoverage();
	if (S==0) return 0;	

	count_t M_=pro.getcount(0);
	count_t m_=pro.getcount(1);

	count_t E_=S-M_-m_;

	float_t e_;

	if (minerr>0){
		e_= float_t(E_)/float_t(S);
		a.error=3.*e_/2.;
		e_=float_t(E_+m_)/float_t(S);
		a.null_error=e_;
		if (a.error<minerr) a.error=minerr;
		if (a.null_error<minerr) a.null_error=minerr;
	} else {
		e_= float_t(E_)/float_t(S);
		a.error=3.*e_/2.;
		e_=float_t(E_+m_)/float_t(S);
		a.null_error=e_;
	};

	a.major=pro.getindex(0);
	a.minor=pro.getindex(1);

//	std::cout << M_ << " (8) " << m_ << std::endl;
	count_t MM_, Mm_, mm_;
	M_=0;m_=0;MM_=0;Mm_=0;mm_=0;
	while (it!=end){
		if (!it->masked){
			float_t T=float_t(count(*it));
			if (T>MIN){
				N_+=1.;
//				std::cout << (*it).base[a.major] << " -- " << (*it).base[a.minor] << std::endl;
//				std::cout << (*it).base[2] << " -- " << (*it).base[3] << std::endl;
				if (float_t( (*it).base[a.major])/T>0.9) {M_+=2.; MM_+=1.;}
				else if (float_t( (*it).base[a.minor])/T>0.9) {m_+=2.; mm_+=1.;}
				else {M_+=1.; m_+=1.; Mm_+=1.;}
			}
		}
		++it;
	}
//	std::cout << MM_ << " |" << Mm_ << "| " << mm_ << std::endl;
	if (N_==0) return 0;	
	float_t P_=MM_, Q_=mm_, H_=Mm_;
	a.N=N_;
	N_=P_+Q_+H_;
	P_/=N_;
	Q_/=N_;
	H_/=N_;
	a.freq=P_+H_/2;
	a.MM=P_;
	a.Mm=H_;
	a.mm=Q_;
	a.efc=efc(pro, MIN);
	return 0;
}


/*TODO Implement!! Uses a sacant method to maximize the likelihood equations. This method has not been implemented  */

/*count_t maximizesecant (profile &pro, allele_stat &a, count_t const &MIN, float_t const &maxgof, count_t const &maxpitch){

	float_t dll1[3], dll2[3], R;
	dll(pro, a, MIN, dll1);
	while (R>0){
		dll(pro, a, MIN, dll2);
		a.MM+=0;
		a.Mm+=0;
		a.mm+=0;
		R=0;
	};
};*/

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximizegrid (profile &pro, allele_stat &a, count_t const &MIN, float_t const &maxgof, count_t const &maxpitch){

	count_t N_=a.N;
	count_t P_=a.MM*N_;
	count_t H_=a.Mm*N_;
	count_t Q_=N_-P_-H_;
	count_t iP, iQ, iH;
	float_t PtQ, PtH, QtP, QtH, HtP, HtQ, maxll_;
	count_t it=0;
	count_t excluded=0;

	bool running=true;

	a.MM=float_t(P_)/float_t(N_);
	a.Mm=float_t(H_)/float_t(N_);
	a.mm=float_t(Q_)/float_t(N_);

	maxll_=loglikelihood(pro, a, MIN);

	while (running){
		it+=1;
		iP=0;
		iQ=0;
		iH=0;
//		std::cout << a.MM << ", " << a.Mm << ", " <<  a.mm <<std::endl;
		if (P_>0){
			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			PtQ=loglikelihood(pro, a, MIN);

			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			PtH=loglikelihood(pro, a, MIN);

			if (PtQ>PtH) { if (PtQ>maxll_) { iP=-1; iQ=1; iH=0; maxll_=PtQ;} }
			else { if (PtH>maxll_){iP=-1;iQ=0;iH=1; maxll_=PtH;} };
		}
		if (Q_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtP=loglikelihood(pro, a, MIN);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtH=loglikelihood(pro, a, MIN);

			if (QtP>QtH) { if (QtP>maxll_) {iP=1;iQ=-1;iH=0;maxll_=QtP;} }
			else { if (QtH>maxll_) {iP=0; iQ=-1; iH=1; maxll_=QtH;} };
		}
		if (H_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			HtP=loglikelihood(pro, a, MIN);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			HtQ=loglikelihood(pro, a, MIN);

			if (HtP>HtQ) { if (HtP>maxll_) {iP=1; iQ=0; iH=-1; maxll_=HtP;} }
			else { if (HtQ>maxll_) {	iP=0; iQ=1; iH=-1; maxll_=HtQ;} }
		}
		if (iP==0 && iQ==0 && iH==0) {running=false;}
		else { P_+=iP; Q_+=iQ; H_+=iH;}
	}

	a.MM=float_t(P_)/float_t(N_);
	a.Mm=float_t(H_)/float_t(N_);
	a.mm=float_t(Q_)/float_t(N_);

	a.freq=float_t(P_)/float_t(N_)+float_t(H_)/float_t(N_)/2.;

	if (a.freq!=0. && a.freq!=1.) a.f=1-float_t(H_)/float_t(N_)/(2.*(1.-a.freq)*a.freq);
	else a.f=0;
	a.gof=gof(pro, a, MIN, maxgof);

	excluded=pro.maskedcount();
	
	if ( abs(a.gof)>maxgof) {
		if (excluded==maxpitch) return excluded;
		return maximizegrid(pro, a, MIN, maxgof, maxpitch);
	};

	return excluded;
}
