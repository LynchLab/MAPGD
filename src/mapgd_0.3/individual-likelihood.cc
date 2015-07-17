#include "individual-likelihood.h"
#include <list>
/* Written by Matthew Ackerman.*/

/* Currently likelihood calculations are all performed by a set of 'models' that we give to the "multinomial" class via the 
   multinomial::set method. These models should be able to look at the allele_stat structure, which contains information about
   the error rate at the locus and the identity of the major and minor allele, and return a set of four log probabilities of 
   observing each particular nucleotide in a given call.
*/

//I'm overloading the += operator so that I can write A+=B, because I'm lazy.

models::models(void){
	lnMM_=new lnmultinomial(4);
	lnMm_=new lnmultinomial(4);
	lnmm_=new lnmultinomial(4);
	lnF_=new lnmultinomial(4);
}

models::~models(void){
	delete lnMM_;
	delete lnMm_;
	delete lnmm_;
	delete lnF_;
}

void MMmodelP(const allele_stat &a, float_t *l){ 	//The model used for the goodness of fit test assuming the genotype is 
						 	//[M]ajor [M]ajor.
	float_t H=(1-a.error);
	float_t e3=log(1-H);
	if (H==1) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
	if (H==0) l[a.major]=-FLT_MAX;
};

void mmmodelP(const allele_stat &a, float_t *l){	//Dito, assuming [m]ajor [m]inor.
	float_t H=(a.error/3);
	float_t e3=log(1-H);
	if (H==1) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
	if (H==0) l[a.major]=-FLT_MAX;
};

void MmmodelP(const allele_stat &a, float_t *l){ 	//[M]ajor [m]inor.
	float_t H=(0.5*(1.-a.error)+0.5*a.error/3.);
	float_t e3=log(1-H);
	if (H==1) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(H);
	if (H==0) e3=-FLT_MAX;
};

void MMmodel(const allele_stat &a, float_t *l){		//These are the full probabilities used for fitting the maximum
	float_t e3=log(a.error/3.);			//likelihood model.
	if (a.error==0) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.major]=log(1.-a.error);
	if (a.error==1.) l[a.major]=-FLT_MAX;
};

void mmmodel(const allele_stat &a, float_t *l){		//Dito.
	float_t e3=log(a.error/3.);
	if (a.error==0.) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	l[a.minor]=log(1.-a.error);
	if (a.error==1.) l[a.minor]=-FLT_MAX;
};

void Mmmodel(const allele_stat &a, float_t *l){		//Dito.
	float_t e3=log(a.error/3.);
	if (a.error==0.) e3=-FLT_MAX;
	l[0]=e3; l[1]=e3; l[2]=e3; l[3]=e3;
	float_t H=log(0.5*(1.-a.error)+0.5*a.error/3.);
	l[a.major]=H;
	l[a.minor]=H;
	if (a.error==1.){
		l[a.major]=-FLT_MAX;
		l[a.minor]=-FLT_MAX;
	};
};

/*! \Breif, a function that calculats the log likelihood of a set of observations. */
float_t models::loglikelihood(const Locus &site, const allele_stat &p, const count_t &MIN){

	float_t sumll=0;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	//Tells us when to stop iterating, so we don't generate a seg. fault.

	lnMM_->set(&MMmodel, p);	//Lets initialize the multinomial distributions that will tell us the probability of  
	lnMm_->set(&Mmmodel, p);	//observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm_->set(&mmmodel, p);	//

	float_t logMM=log(p.MM);	//The frequency of the MM genotype in the sample.
	float_t logMm=log(p.Mm);	// Dito Mm.
	float_t logmm=log(p.mm);	// Diot mm.

	float_t E0, E1, E2;		//These are some variables to breifly store a portion of our likelihood calculation

	count_t T;

	while(it!=end){
		if (!it->masked){
			T=( it->base[0]+it->base[1]+it->base[2]+it->base[3] );
			if ( T>=MIN ) {

				E0=logMM+lnMM_->lnprob(it->base);
				E1=logMm+lnMm_->lnprob(it->base);
				E2=logmm+lnmm_->lnprob(it->base);

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
	if(isnan(sumll) ) sumll=-FLT_MAX;
	return sumll;
}

/*! \Breif, a function that calculats the log likelihood of a set of observations. */
float_t models::genotypelikelihood(quartet_t const &quartet, const allele_stat &population, const count_t &MIN){

	lnMM_->set(&MMmodel, population);	//Lets initialize the multinomial distributions that will tell us the probability of  

	lnMm_->set(&Mmmodel, population);  //observing a particular quartet given that the individual has the MM, Mm or mm genotype.
	lnmm_->set(&mmmodel, population);  //

	float_t logMM=log(population.MM); //The frequency of the MM genotype in the sample.
	float_t logMm=log(population.Mm); // Dito Mm.
	float_t logmm=log(population.mm); // Diot mm.

	float_t E[3];

	E[0]=logMM+lnMM_->lnprob(quartet.base);
	E[1]=logMm+lnMm_->lnprob(quartet.base);
	E[2]=logmm+lnmm_->lnprob(quartet.base);

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
	return E[0];
}

/*@Breif, Used for GOF, same basic idea as above. I know, if I have to write it twice I'm doing something wrong.*/
float_t models::lnP(const count_t *base, const allele_stat &p){

	lnMM_->set(&MMmodelP, p);
	lnMm_->set(&MmmodelP, p);
	lnmm_->set(&mmmodelP, p);
	
	float_t logMM=log(p.MM);
	float_t logMm=log(p.Mm);
	float_t logmm=log(p.mm);

	float_t E0, E1, E2, tE;

	E0=logMM+lnMM_->lnprob(base);
	E1=logMm+lnMm_->lnprob(base);
	E2=logmm+lnmm_->lnprob(base);

	if(E0>E2) {tE=E2; E2=E0; E0=tE;}
	if(E1>E2) {tE=E2; E2=E1; E1=tE;}

	return log(exp(E0-E2)+exp(E1-E2)+1.)+E2;
};

/*Used for GOF, this calculate the expectation ad variance. The majority of the range of x may contribute essentially nothing 
  to this calculation, so I may look into some better bounds to speed up these calculations.*/

void EVofP (count_t N_, const allele_stat &a, models &model, float_t &E, float_t &V){
	count_t l[4];
	float_t tP;
	E=0; V=0;

	//E=Sum[(Binomial[N, x]^2 p^(2 x) (1 - p)^(2 (N - x)))/N, {x, 0, N}]
	//E=((p - 1)^(2 N) H2F1[-N, -N, 1, p^2/(p-1)^2])/N
	//V=(E^2-2*(1-p)^(2*N)*E*H2F1[-N, -N, 1, p^2/(p-1)^2 ]+(1-p)^(3*N)*H2F1PFQ[{-N, -N, -N}, {1, 1}, p^3/(p-1)^3 ] )/N
	
//	float_t A[3]={-N_, -N_, -N_};
//	float_t B[2]={1., 1.};
//	float_t p=
//	std::cout << "---===---\n";
//	E=(pow(p-1, 2*N_)*pFq::_2F1(-N_, -N_, 1., pow(p, 2) /pow(p-1,2) )/N_;
	float_t H=(0.5*(1.-a.error)+0.5*a.error/3.);
//	N_=10;
	std::cout <<  "E " << a.error << ", " << N_ << "\n" <<
		a.MM*pow(1.-a.error-1., 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(1-a.error, 2) /pow( (1.-a.error)-1.,2))
		+a.mm*pow(a.error/3.-1, 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(a.error/3., 2) /pow(a.error/3.-1.,2))
		+a.Mm*pow(H-1., 2.*N_)*pFq::_2F1(-float_t(N_), -float_t(N_), 1., pow(H, 2) /pow(H-1.,2))
		 << '\n';
//	std::cout << (pow(E, 2)-pow(2*(1-p), 2*N)*E*pFq::_2F1( -N_, -N_, 1., pow(p, 2)/pow(p-1, 2) )+pow(1-p, 3*N)*pFq::_3F2(A, B, pow(p, 3)/pow(p-1, 3) ) )/N_ << '\n';
//	E=0;
	for (size_t x=0; x<N_+1; ++x){
		memset(l, 0, sizeof(count_t)*4);
		l[a.major]=x;
		l[a.minor]=N_-x;
		tP=model.lnP(l, a);
		E+=exp(tP)*exp(tP);
		V+=exp(tP)*pow(tP, 2);
	}
	V-=pow(E, 2);
	std::cout << "O " << E << "\n\n";
//	std::cout << V << '\n';
}

//THE goodness of fit calcualtion.
float_t gof (Locus &site, const allele_stat &a, models &model, std::vector <float_t> &gofs, const count_t &MIN, const float_t &site_maxgof){
	float_t Num=0., Den=0., E, V, O, thisgof,clone_mingof=FLT_MAX;
	std::vector <float_t>::iterator gof=gofs.begin();

	count_t M_, N_, l[4];

	std::vector <quartet_t>::iterator it=site.sample.begin(); 
	std::vector <quartet_t>::iterator end=site.sample.end(); 

	quartet_t * maxgof_ptr; 

	while (it!=end){
		if (!it->masked){
			N_=float_t(count(*it));
			if ( N_>MIN ){
				memset(l, 0, sizeof(count_t)*4);
				M_=(*it).base[a.major];
				l[a.major]=M_;
				l[a.minor]=N_-M_;
				//
				EVofP(N_, a, model, E, V);
				//
				O=model.lnP(l, a);
				Num+=O-E;
				Den+=V;
				if (V>0) (*gof)=( (O-E)/sqrt(V) );
				if((*gof)<clone_mingof){
					clone_mingof=(*gof);
					maxgof_ptr=&(*it);
				};
			}
		}
		++it;
		++gof;
	}
	if (Den<=0) return 0.;
	thisgof=Num/sqrt(Den);
	if(abs(thisgof)>site_maxgof) site.mask(maxgof_ptr);
	return thisgof;
}

/*calculates the 'effective number of chromosomes' in a sample. This number isn't used for anything right now.*/

float_t efc (const Locus &site, const count_t &MIN){
	std::vector <quartet_t>::const_iterator it=site.sample.begin(); 
	std::vector <quartet_t>::const_iterator end=site.sample.end(); 

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

//Intilaizes the parameters of the ... returns 0 on succesful excecution, returns 1 if there is am...
count_t init_params(Locus &site, allele_stat &a, const count_t &MIN, const float_t &minerr, const count_t &M){

	std::vector <quartet_t>::iterator it=site.sample.begin(); 
	std::vector <quartet_t>::iterator end=site.sample.end(); 

	site.sort();

	count_t N_=0;

	count_t S=site.getcoverage();

	if (S==0) {
		a.major=site.getindex(0);
		a.minor=site.getindex(1);
		return 0;
	}

	count_t M_=site.getcount(0);
	count_t m_=site.getcount(1);

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

	a.major=site.getindex(0);
	a.minor=site.getindex(1);
	a.e1=site.getindex(2);
	a.e2=site.getindex(3);

	count_t MM_, Mm_, mm_;
	M_=0;m_=0;MM_=0;Mm_=0;mm_=0;
	while (it!=end){
		if (!it->masked){
			float_t T=float_t(count(*it));
			if (T>MIN){
				N_+=1.;
				if (float_t( (*it).base[a.major])/T>0.9) {M_+=2.; MM_+=1.;}
				else if (float_t( (*it).base[a.minor])/T>0.9) {m_+=2.; mm_+=1.;}
				else {M_+=1.; m_+=1.; Mm_+=1.;}
			}
			else { it->masked=true;};
		}
		++it;
	}
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
	a.efc=efc(site, MIN);
	return 1;
}


count_t maximize_newton (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const count_t &MIN, const float_t &maxgof, const count_t &maxpitch){

	float_t J[3][3];		//The Jacobian/Hessian.
        float_t iJ[3][3]; 		//Inverse of the Jacobian/Hessian.
	float_t R[3]={100,100,100};	//The 'residual'. Th

        float_t det, lim;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	 

	count_t iter=0;			//counts the number of iterations to let us know if we have a failure to converge.

        while ( ( (fabs(R[0])+fabs(R[1])+fabs(R[2]) )>0.00001 || isnan(R[0]) || isnan(R[1]) || isnan(R[2]) ) && iter<100){

		++iter;
 
		memset(J[0], 0, sizeof(float_t)*3);
		memset(J[1], 0, sizeof(float_t)*3);
		memset(J[2], 0, sizeof(float_t)*3);
		memset(R, 0, sizeof(float_t)*3);

		it=site.sample.begin();	 

		while (it!=end ){

			J[0][0]+=J00(*it, a); J[0][1]+=J01(*it, a); J[0][2]+=J02(*it, a);
			J[1][0]+=J10(*it, a); J[1][1]+=J11(*it, a); J[1][2]+=J12(*it, a);
			J[2][0]+=J20(*it, a); J[2][1]+=J21(*it, a); J[2][2]+=J22(*it, a);

                        R[0]+=H0(*it, a);
                        R[1]+=H1(*it, a);
                        R[2]+=H2(*it, a);

                        ++it;
                };
/*		Check my matrix algebra.
		
		std::cout << "==\n";

		J[0][0]=0; J[0][1]=3; J[0][2]=6;
		J[1][0]=1; J[1][1]=4; J[1][2]=7;
		J[2][0]=1; J[2][1]=5; J[2][2]=8;
						*/

		//FIRST ROW 
		iJ[0][0]= (J[1][1]*J[2][2]-J[2][1]*J[1][2]); 
		iJ[0][1]=-(J[0][1]*J[2][2]-J[0][2]*J[2][1]); 
		iJ[0][2]= (J[0][1]*J[1][2]-J[0][2]*J[1][1]);

		//SECOND ROW 
		iJ[1][0]=-(J[1][0]*J[2][2]-J[1][2]*J[2][0]); 
		iJ[1][1]= (J[0][0]*J[2][2]-J[0][2]*J[2][0]); 
		iJ[1][2]=-(J[0][0]*J[1][2]-J[0][2]*J[1][0]);

		//THIRD ROW
                iJ[2][0]= (J[1][0]*J[2][1]-J[1][1]*J[2][0]); 
		iJ[2][1]=-(J[0][0]*J[2][1]-J[0][1]*J[2][0]); 
		iJ[2][2]= (J[0][0]*J[1][1]-J[0][1]*J[1][0]);

		//THE DETERMINENT
		det=J[0][0]*iJ[0][0]+J[0][1]*iJ[1][0]+J[0][2]*iJ[2][0];

		iJ[0][0]/=det; iJ[0][1]/=det; iJ[0][2]/=det;
		iJ[1][0]/=det; iJ[1][1]/=det; iJ[1][2]/=det;
		iJ[2][0]/=det; iJ[2][1]/=det; iJ[2][2]/=det;

                R[0]=(R[0]*iJ[0][0]+R[1]*iJ[0][1]+R[2]*iJ[0][2]);
                R[1]=(R[0]*iJ[1][0]+R[1]*iJ[1][1]+R[2]*iJ[1][2]);
                R[2]=(R[0]*iJ[2][0]+R[1]*iJ[2][1]+R[2]*iJ[2][2]);
	
//	       	std::cout << "Pi=" << a.freq-R[0] << ", Epsilon=" << a.error-R[1]  << ", F=" << a.f-R[2] << ", R=" << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ": " << model.loglikelihood(site, a, MIN) << std::endl;

		//BONDS CHECKING.
		if (a.freq-R[0]>1.0) a.freq=1-(1-a.freq)/2.;
                else if (a.freq-R[0]<0) a.freq/=2.0;
		else a.freq-=R[0];

		if (a.error-R[1]>0.25) a.error=0.25-(0.25-a.error)/2.;
                else if (a.error-R[1]<0) a.error/=2.0;
		else a.error-=R[1];

		lim=-pow( (1-a.freq), 2)/(a.freq*(1-a.freq) );
		if (a.f-R[2]>1.0) a.f=1-(1-a.f)/2.;
                else if (a.f-R[2]<lim) a.f=lim+(lim-a.f)/2.;
		else a.f-=R[2];
//	       	std::cout << "Pi=" << a.freq << ", Epsilon=" << a.error  << ", F=" << a.f << ", R=" << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ": " << model.loglikelihood(site, a, MIN) << std::endl;
        };

        a.MM=pow(a.freq,2)+a.f*(a.freq*(1-a.freq) );
        a.Mm=2*(1-a.f)*(a.freq*(1-a.freq) );
        a.mm=pow(1-a.freq,2)+a.f*(a.freq*(1-a.freq) );
	a.coverage=site.getcoverage();

	count_t excluded=0;

	excluded=site.maskedcount();
	std::vector <float_t> temp_gofs(site.sample.size());
	a.gof=gof(site, a, model, temp_gofs, MIN, maxgof);
	if (iter==100) {
		std::cerr << "Failure to maximize " << a << "\n";
		init_params(site, a, MIN, 0, 0);
		return maximize_grid(site, a, model, gofs, MIN, maxgof, maxpitch);
	}
	if ( abs(a.gof)>maxgof) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
		return maximize_newton(site, a, model, gofs, MIN, maxgof, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return excluded;
}

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximize_grid (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const count_t &MIN, const float_t &maxgof, const count_t &maxpitch){

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
	
	maxll_=model.loglikelihood(site, a, MIN);

	while (running){
		it+=1;
		iP=0;
		iQ=0;
		iH=0;
		if (P_>0){
			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			PtQ=model.loglikelihood(site, a, MIN);

			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			PtH=model.loglikelihood(site, a, MIN);

			if (PtQ>PtH) { if (PtQ>maxll_) { iP=-1; iQ=1; iH=0; maxll_=PtQ;} }
			else { if (PtH>maxll_){iP=-1;iQ=0;iH=1; maxll_=PtH;} };
		}
		if (Q_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtP=model.loglikelihood(site, a, MIN);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtH=model.loglikelihood(site, a, MIN);

			if (QtP>QtH) { if (QtP>maxll_) {iP=1;iQ=-1;iH=0; maxll_=QtP;} }
			else { if (QtH>maxll_) {iP=0; iQ=-1; iH=1; maxll_=QtH;} };
		}
		if (H_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			HtP=model.loglikelihood(site, a, MIN);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			HtQ=model.loglikelihood(site, a, MIN);
			
			if (HtP>HtQ) { if (HtP>maxll_) {iP=1; iQ=0; iH=-1; maxll_=HtP;} }
			else { if (HtQ>maxll_) {	iP=0; iQ=1; iH=-1; maxll_=HtQ;} }
		}
		if (iP==0 && iQ==0 && iH==0) {running=false;}
		else { P_+=iP; Q_+=iQ; H_+=iH;}
	}

	a.coverage=site.getcoverage();

	a.MM=float_t(P_)/float_t(N_);
	a.Mm=float_t(H_)/float_t(N_);
	a.mm=float_t(Q_)/float_t(N_);

	a.freq=float_t(P_)/float_t(N_)+float_t(H_)/float_t(N_)/2.;

	if (a.freq!=0. && a.freq!=1.) a.f=1-float_t(H_)/float_t(N_)/(2.*(1.-a.freq)*a.freq);
	else a.f=0;

	std::vector <float_t> temp_gofs(site.sample.size());
	a.gof=gof(site, a, model, temp_gofs, MIN, maxgof);

	excluded=site.maskedcount();
	
	if ( abs(a.gof)>maxgof) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
		return maximize_grid(site, a, model, gofs, MIN, maxgof, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return excluded;
}

count_t maximize_analytical (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const count_t &MIN, const float_t &maxgof, const count_t &maxpitch){

	a.coverage=site.getcoverage();

	a.MM=1.0;
	a.Mm=0.0;
	a.mm=0.0;

	a.freq=1.0;

	//These are both undefined...
	a.f=0;
	a.gof=0;		

	return 0;
}
