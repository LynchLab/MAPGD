#include "individual-likelihood.h"
#include <list>
/* Written by Matthew Ackerman.*/

/*calculates the 'effective number of chromosomes' in a sample. This number isn't used for anything right now.*/

float_t efc (const Locus &site){
	std::vector <quartet_t>::const_iterator it=site.sample.begin(); 
	std::vector <quartet_t>::const_iterator end=site.sample.end(); 

	float_t ec=0;

	while (it!=end){
		if (!it->masked){
			ec+=2.-2.*pow(0.5, count(*it) );
		}
		++it;
	}
	return ec;
};

//Intilaizes the parameters of the ... returns 0 on succesful excecution, returns 1 if there is am...
count_t init_params(Locus &site, allele_stat &a, const float_t &minerr){

	a.coverage=site.getcoverage();

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
	if (m_==site.getcount(2) ){
		if (m_==site.getcount(3) ) a.minor=site.getindex(1+rand()%3);
		else a.minor=site.getindex(1+rand()%2);
	};

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
			N_+=1.;
			if (float_t( (*it).base[a.major])/T>0.9) {M_+=2.; MM_+=1.;}
			else if (float_t( (*it).base[a.minor])/T>0.9) {m_+=2.; mm_+=1.;}
			else {M_+=1.; m_+=1.; Mm_+=1.;}
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
	return 1;
}


count_t maximize_newton (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const float_t &maxgof, const count_t &maxpitch){

	float_t J[3][3];		//The Jacobian/Hessian.
        float_t iJ[3][3]; 		//Inverse of the Jacobian/Hessian.
	float_t R[3]={100,100,100};	//The 'residual'. Th

        float_t det;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	 

	count_t iter=0;			//counts the number of iterations to let us know if we have a failure to converge.

	if(a.error==0) a.error=0.005;

        while ( ( (fabs(R[0])+fabs(R[1])+fabs(R[2]) )>0.00001 || isnan(R[0]) || isnan(R[1]) || isnan(R[2]) ) && iter<200){

		++iter;
 
		memset(J[0], 0, sizeof(float_t)*3);
		memset(J[1], 0, sizeof(float_t)*3);
		memset(J[2], 0, sizeof(float_t)*3);
		memset(R, 0, sizeof(float_t)*3);

		it=site.sample.begin();	 

		float_t sumlnL=0;

		while (it!=end ){

			J[0][0]+=J00(*it, a); J[0][1]+=J01(*it, a); J[0][2]+=J02(*it, a);
			J[1][0]+=J10(*it, a); J[1][1]+=J11(*it, a); J[1][2]+=J12(*it, a);
			J[2][0]+=J20(*it, a); J[2][1]+=J21(*it, a); J[2][2]+=J22(*it, a);

                        R[0]+=H0(*it, a);
                        R[1]+=H1(*it, a);
                        R[2]+=H2(*it, a);
			
			sumlnL+=lnL_NR(*it, a);

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

//		std::cout << R[0] << "*" << iJ[0][0] << std::endl;

                R[0]=(R[0]*iJ[0][0]+R[1]*iJ[0][1]+R[2]*iJ[0][2]);
                R[1]=(R[0]*iJ[1][0]+R[1]*iJ[1][1]+R[2]*iJ[1][2]);
                R[2]=(R[0]*iJ[2][0]+R[1]*iJ[2][1]+R[2]*iJ[2][2]);
	
	       	std::cout << iter <<" Pi= " << a.freq << ", Epsilon=" << a.error  << ", F=" << a.f << ", R=" << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ": " << sumlnL << " : " << sizeof(float_t)  << " : "<< model.loglikelihood(site, a) << std::endl;

		//BOUNDS CHECKS.
		/*
		if (a.freq-R[0]>1.0) a.freq=1-(1-a.freq)/2.;
                else if (a.freq-R[0]<0) a.freq/=2.0;
		else a.freq-=R[0];

		if (a.error-R[1]>0.25) a.error=0.25-(0.25-a.error)/2.;
                else if (a.error-R[1]<0) a.error/=2.0;
		else a.error-=R[1];

		lim=-pow( (1-a.freq), 2)/(a.freq*(1-a.freq) );
		if (a.f-R[2]>1.0) a.f=1-(1-a.f)/2.;
                else if (a.f-R[2]<lim) a.f=lim+(lim-a.f)/2.;
		else a.f-=R[2];*/
		a.freq-=R[0];
		a.error-=R[1];
		a.f-=R[2];
        };

        a.MM=pow(a.freq,2)+a.f*(a.freq*(1-a.freq) );
        a.Mm=2*(1-a.f)*(a.freq*(1-a.freq) );
        a.mm=pow(1-a.freq,2)+a.f*(a.freq*(1-a.freq) );
	a.coverage=site.getcoverage();

	count_t excluded=0;

	excluded=site.maskedcount();
	std::vector <float_t> temp_gofs(site.sample.size());
	a.gof=model.goodness_of_fit(site, a, temp_gofs, maxgof);
	if (iter==200) {
		std::cerr << "Failure to maximize " << iter << " " << a << "\n";
		init_params(site, a, 0);
		return maximize_grid(site, a, model, gofs, maxgof, maxpitch);
	}
	return 0;
	if ( a.gof<maxgof) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
		return maximize_newton(site, a, model, gofs, maxgof, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return excluded;
}

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximize_grid (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const float_t &MINGOF, const size_t &maxpitch){

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
	
	maxll_=model.loglikelihood(site, a);

	while (running){
		it+=1;
		iP=0;
		iQ=0;
		iH=0;
		if (P_>0){
			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			PtQ=model.loglikelihood(site, a);

			a.MM=float_t(P_-1)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			PtH=model.loglikelihood(site, a);

			if (PtQ>PtH) { if (PtQ>maxll_) { iP=-1; iQ=1; iH=0; maxll_=PtQ;} }
			else { if (PtH>maxll_){iP=-1;iQ=0;iH=1; maxll_=PtH;} };
		}
		if (Q_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtP=model.loglikelihood(site, a);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_+1)/float_t(N_);
			a.mm=float_t(Q_-1)/float_t(N_);

			QtH=model.loglikelihood(site, a);

			if (QtP>QtH) { if (QtP>maxll_) {iP=1;iQ=-1;iH=0; maxll_=QtP;} }
			else { if (QtH>maxll_) {iP=0; iQ=-1; iH=1; maxll_=QtH;} };
		}
		if (H_>0){
			a.MM=float_t(P_+1)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_)/float_t(N_);

			HtP=model.loglikelihood(site, a);

			a.MM=float_t(P_)/float_t(N_);
			a.Mm=float_t(H_-1)/float_t(N_);
			a.mm=float_t(Q_+1)/float_t(N_);

			HtQ=model.loglikelihood(site, a);
			
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
	a.gof=model.goodness_of_fit(site, a, temp_gofs, MINGOF);

	excluded=site.maskedcount();
	a.efc=efc(site);
	
	if ( a.gof<MINGOF) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
		a.N-=1;
		return maximize_grid(site, a, model, gofs, MINGOF, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return 	excluded;
}



count_t maximize_analytical (Locus &site, allele_stat &a, models &model, std::vector <float_t> &gofs, const float_t &maxgof, const size_t &maxpitch){

	a.MM=1.0;
	a.Mm=0.0;
	a.mm=0.0;

	a.freq=1.0;

	//These are both undefined...
	a.f=0;
	a.gof=0;		

	return 0;
}
