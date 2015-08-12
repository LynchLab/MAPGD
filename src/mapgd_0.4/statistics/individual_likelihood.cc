#include "individual_likelihood.h"
#include <list>
/* Written by Matthew Ackerman.*/

/*calculates the 'effective number of chromosomes' in a sample. This number isn't used for anything right now.*/

real_t efc (const locus &site){
	std::vector <quartet_t>::const_iterator it=site.begin(); 
	std::vector <quartet_t>::const_iterator end=site.end(); 

	real_t ec=0;
	count_t N_;

	while (it!=end){
		if (!it->masked){
			ec+=2.-2.*pow(0.5, N_);
		}
		++it;
	}
	return ec;
};

//Intilaizes the parameters of the ... returns 1 on succesful excecution, returns 0 if there are no reads...
count_t init_params(locus &site, allele &a, const real_t &minerr){

	std::vector <quartet_t>::iterator it=site.begin(); 
	std::vector <quartet_t>::iterator end=site.end(); 

	site.sort();

	count_t N_=0;

	count_t S=site.get_coverage();

	if (S==0) {
		a.major=site.get_index(0);
		a.minor=site.get_index(1);
		a.error=0;
		return 0;
	}

	count_t M_=site.get_count(0);
	count_t m_=site.get_count(1);

	count_t E_=S-M_-m_;

	real_t e_;

	if (minerr>0){
		e_= real_t(E_)/real_t(S);
		a.error=3.*e_/2.;
		e_=real_t(E_+m_)/real_t(S);
		if (a.error<minerr) a.error=minerr;
	} else {
		e_= real_t(E_)/real_t(S);
		a.error=3.*e_/2.;
		e_=real_t(E_+m_)/real_t(S);
		if (a.error<minerr) a.error=minerr;
	};

	a.base[0]=site.get_index(0);
	a.base[1]=site.get_index(1);
	a.base[2]=site.get_index(2);
	a.base[3]=site.get_index(3);

	count_t MM_, Mm_, mm_;
	M_=0;m_=0;MM_=0;Mm_=0;mm_=0;
	while (it!=end){
		if (!it->masked){
			real_t T=real_t(count(*it));
			N_+=1.;
			if (real_t( (*it).base[a.major])/T>0.9) {M_+=2.; MM_+=1.;}
			else if (real_t( (*it).base[a.minor])/T>0.9) {m_+=2.; mm_+=1.;}
			else {M_+=1.; m_+=1.; Mm_+=1.;}
		}
		++it;
	}
	if (N_==0) return 0;	
	real_t P_=MM_, Q_=mm_, H_=Mm_;
	N_=P_+Q_+H_;
	P_/=N_;
	Q_/=N_;
	H_/=N_;
	a.frequency[0]=P_;
	a.frequency[1]=H_;
	a.frequency[2]=Q_;
	return 1;
}


count_t maximize_newton (locus &site, allele &a, models &model, std::vector <real_t> &gofs, const real_t &maxgof, const count_t &maxpitch){

	real_t J[3][3];		//The Jacobian/Hessian.
        real_t iJ[3][3]; 		//Inverse of the Jacobian/Hessian.
	real_t R[3]={100,100,100};	//The 'residual'. Th

        real_t det, lim;

	std::vector <quartet_t>::const_iterator it=site.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.end();	 

	count_t iter=0;			//counts the number of iterations to let us know if we have a failure to converge.

	if(a.error==0) a.error=0.005;

	return 0;
//        while ( ( (fabs(R[0])+fabs(R[1])+fabs(R[2]) )>0.00001 || isnan(R[0]) || isnan(R[1]) || isnan(R[2]) ) && iter<200){
/*
		++iter;
 
		memset(J[0], 0, sizeof(real_t)*3);
		memset(J[1], 0, sizeof(real_t)*3);
		memset(J[2], 0, sizeof(real_t)*3);
		memset(R, 0, sizeof(real_t)*3);

		it=site.begin();	 

		real_t sumlnL=0;

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
*//*		Check my matrix algebra.
		
		std::cout << "==\n";

		J[0][0]=0; J[0][1]=3; J[0][2]=6;
		J[1][0]=1; J[1][1]=4; J[1][2]=7;
		J[2][0]=1; J[2][1]=5; J[2][2]=8;
						*/

/*		//FIRST ROW 
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
	
	       	std::cout << iter <<" Pi= " << a.freq << ", Epsilon=" << a.error  << ", F=" << a.f << ", R=" << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ": " << sumlnL << " : " << sizeof(real_t)  << " : "<< model.loglikelihood(site, a) << std::endl;

		//BOUNDS CHECKS.
*/		/*
		if (a.freq-R[0]>1.0) a.freq=1-(1-a.freq)/2.;
                else if (a.freq-R[0]<0) a.freq/=2.0;
		else a.freq-=R[0];

		if (a.error-R[1]>0.25) a.error=0.25-(0.25-a.error)/2.;
                else if (a.error-R[1]<0) a.error/=2.0;
		else a.error-=R[1];

		lim=-pow( (1-a.freq), 2)/(a.freq*(1-a.freq) );
		if (a.f-R[2]>1.0) a.f=1-(1-a.f)/2.;
                else if (a.f-R[2]<lim) a.f=lim+(lim-a.f)/2.;
		else a.f-=R[2];*/ /*
		a.freq-=R[0];
		a.error-=R[1];
		a.f-=R[2];
*/ /*     };

        a.MM=pow(a.freq,2)+a.f*(a.freq*(1-a.freq) );
        a.Mm=2*(1-a.f)*(a.freq*(1-a.freq) );
        a.mm=pow(1-a.freq,2)+a.f*(a.freq*(1-a.freq) );
	a.coverage=site.get_coverage();

	count_t excluded=0;

	excluded=site.maskedcount();
	std::vector <real_t> temp_gofs(site.quartet_size());
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
	return excluded;*/
}

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximize_grid (locus &site, allele &a, models &model, std::vector <real_t> &gofs, const real_t &MINGOF, const size_t &maxpitch){

	count_t N_=site.get_coverage();
	count_t P_=a.MM*N_;
	count_t H_=a.Mm*N_;
	count_t Q_=N_-P_-H_;
	count_t iP, iQ, iH;
	real_t PtQ, PtH, QtP, QtH, HtP, HtQ, maxll_;
	count_t it=0;
	count_t excluded=0;

	bool running=true;

	a.MM=real_t(P_)/real_t(N_);
	a.Mm=real_t(H_)/real_t(N_);
	a.mm=real_t(Q_)/real_t(N_);
	
	maxll_=model.loglikelihood(site, a);

	while (running){
		it+=1;
		iP=0;
		iQ=0;
		iH=0;
		if (P_>0){
			a.MM=real_t(P_-1)/real_t(N_);
			a.Mm=real_t(H_)/real_t(N_);
			a.mm=real_t(Q_+1)/real_t(N_);

			PtQ=model.loglikelihood(site, a);

			a.MM=real_t(P_-1)/real_t(N_);
			a.Mm=real_t(H_+1)/real_t(N_);
			a.mm=real_t(Q_)/real_t(N_);

			PtH=model.loglikelihood(site, a);

			if (PtQ>PtH) { if (PtQ>maxll_) { iP=-1; iQ=1; iH=0; maxll_=PtQ;} }
			else { if (PtH>maxll_){iP=-1;iQ=0;iH=1; maxll_=PtH;} };
		}
		if (Q_>0){
			a.MM=real_t(P_+1)/real_t(N_);
			a.Mm=real_t(H_)/real_t(N_);
			a.mm=real_t(Q_-1)/real_t(N_);

			QtP=model.loglikelihood(site, a);

			a.MM=real_t(P_)/real_t(N_);
			a.Mm=real_t(H_+1)/real_t(N_);
			a.mm=real_t(Q_-1)/real_t(N_);

			QtH=model.loglikelihood(site, a);

			if (QtP>QtH) { if (QtP>maxll_) {iP=1;iQ=-1;iH=0; maxll_=QtP;} }
			else { if (QtH>maxll_) {iP=0; iQ=-1; iH=1; maxll_=QtH;} };
		}
		if (H_>0){
			a.MM=real_t(P_+1)/real_t(N_);
			a.Mm=real_t(H_-1)/real_t(N_);
			a.mm=real_t(Q_)/real_t(N_);

			HtP=model.loglikelihood(site, a);

			a.MM=real_t(P_)/real_t(N_);
			a.Mm=real_t(H_-1)/real_t(N_);
			a.mm=real_t(Q_+1)/real_t(N_);

			HtQ=model.loglikelihood(site, a);
			
			if (HtP>HtQ) { if (HtP>maxll_) {iP=1; iQ=0; iH=-1; maxll_=HtP;} }
			else { if (HtQ>maxll_) {	iP=0; iQ=1; iH=-1; maxll_=HtQ;} }
		}
		if (iP==0 && iQ==0 && iH==0) {running=false;}
		else { P_+=iP; Q_+=iQ; H_+=iH;}
	}

	a.MM=real_t(P_)/real_t(N_);
	a.Mm=real_t(H_)/real_t(N_);
	a.mm=real_t(Q_)/real_t(N_);

	std::vector <real_t> temp_gofs(site.quartet_size());
//	a.gof=model.goodness_of_fit(site, a, temp_gofs, MINGOF);

	excluded=site.maskedcount();
	
/*	if ( gof<MINGOF) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
		return maximize_grid(site, a, model, gofs, MINGOF, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
*/	return excluded;
}

count_t maximize_analytical (locus &site, allele &a, models &model, std::vector <real_t> &gofs, const real_t &maxgof, const size_t &maxpitch){


	a.frequency[0]=1.0;
	a.frequency[1]=0.0;
	a.frequency[2]=0.0;

	//These are both undefined...

	return 0;
}
