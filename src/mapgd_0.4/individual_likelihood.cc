#include "individual_likelihood.h"
#include <list>

/* Written by Matthew Ackerman.*/

/*calculates the 'effective number of chromosomes' in a sample. This number isn't used for anything right now.*/

float_t efc (const Locus &site)
{
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
}

//Intilaizes the parameters of the ... returns 0 on succesful excecution, returns 1 if there is am...
count_t init_params(Locus &site, Allele &a, const float_t &minerr)
{
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

	e_= float_t(E_)/float_t(S);
	a.error=3.*e_/2.;
	e_=float_t(E_+m_)/float_t(S);
	a.null_error=e_;
	e_=float_t(E_+M_)/float_t(S);
	a.null_error2=e_;

	if (a.error<minerr) a.error=minerr;
	if (a.null_error<minerr) a.null_error=minerr;
	if (a.null_error2<minerr) a.null_error2=minerr;

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

count_t maximize_restricted_newton (Locus &site, Allele &a, models &model, std::vector <float_t> &gofs, const float_t &maxgof, const size_t &maxpitch)
{
	float_t J[2][2];		//The Jacobian/Hessian.
        float_t iJ[2][2]; 		//Inverse of the Jacobian/Hessian.
	float_t R[2]={100,100};	//The 'residual'. Th

        float_t det;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	 

	count_t iter=0;			//counts the number of iterations to let us know if we have a failure to converge.

	float_t MM=a.MM, Mm=a.Mm, mm=a.mm, p=a.freq, F=a.freq, E;//, q=1-a.freq;
//	a.f=log( (1.+p/(1.-p) )/(a.f+p/(1-p) )-1.);
//	a.f=log( (1.+(1.-p)/p )/(0.2+(1.-p)/p )-1.);
//	MM=pow(p, 2)+p*q*F;
//	Mm=2*p*q*(1-F);
//	mm=pow(q, 2)+p*q*F;

/*
	if (mm==0){
		if (Mm>MM) {mm+=0.01; Mm-=0.01;}
		else {mm+=0.01; MM-=0.01;}
	};
	if (Mm==0){
		if (mm>MM) {Mm+=0.01; mm-=0.01;}
		else {Mm+=0.01; MM-=0.01;}
	};
	if (MM==0){
		if (Mm>mm) {MM+=0.01; Mm-=0.01;}
		else {MM+=0.01; mm-=0.01;}
	};*/

	a.freq=Mm;
	a.f=MM/(MM+mm);


	if (a.freq==1.) a.freq-=0.000001;
	if (a.f==1.) a.f-=0.000001;
	if (a.freq==0.) a.freq+=0.000001;
	if (a.f==0.) a.f+=0.000001;

	a.f=log(1./a.f-1.);
	a.freq=log(1./a.freq-1.);
	a.error=log(0.75/a.error-1.);

	float_t deltalnL=10, maxlnL=DBL_MAX;

        while ( ( (fabs(R[0])+fabs(R[1]) )>0.00001 || std::isnan(R[0]) || std::isnan(R[1]) ) && iter<200 && fabs(deltalnL)>0.0001 ){
#ifdef DEBUG
		std::cerr << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ", " << iter << ", " << fabs(deltalnL) << std::endl;
#endif
		++iter;
 
		memset(J[0], 0, sizeof(float_t)*2);
		memset(J[1], 0, sizeof(float_t)*2);
		memset(R, 0, sizeof(float_t)*2);

		it=site.sample.begin();	 

		float_t sumlnL=0;

		while (it!=end ){

			J[0][0]+=J00R(*it, a); J[0][1]+=J01R(*it, a);
			J[1][0]+=J10R(*it, a); J[1][1]+=J11R(*it, a); 

                        R[0]+=H0R(*it, a);
                        R[1]+=H1R(*it, a);
			
			sumlnL+=lnL_NR(*it, a);

                        ++it;
                };
		deltalnL=(maxlnL-sumlnL)+1;
		maxlnL=sumlnL;

/*		Check my matrix algebra.
		
		J[0][0]=0; J[0][1]=3; J[0][2]=6;
		J[1][0]=1; J[1][1]=4; J[1][2]=7;
		J[2][0]=1; J[2][1]=5; J[2][2]=8;
						*/

		//FIRST ROW 
		iJ[0][0]=J[1][1]; 
		iJ[0][1]=-J[0][1]; 

		//SECOND ROW 
		iJ[1][0]=-J[1][0]; 
		iJ[1][1]= J[0][0]; 

		//THE DETERMINENT
		det=J[0][0]*iJ[0][0]+J[0][1]*iJ[1][0];

		if (fabs(det)<0.00001) 
		{
			det < 0 ? det=-0.00001 : det=0.00001;
		}

		iJ[0][0]/=det; iJ[0][1]/=det;
		iJ[1][0]/=det; iJ[1][1]/=det; 

                R[0]=(R[0]*iJ[0][0]+R[1]*iJ[0][1]);
                R[1]=(R[0]*iJ[1][0]+R[1]*iJ[1][1]);

		p=a.freq;
		F=a.f;

	        a.freq=1./(1.+exp(a.freq) );
	        a.f=1./(1.+exp(a.f) );

	        a.MM=(1.-a.freq)*a.f;
	        a.Mm=a.freq;
	        a.mm=(1.-a.freq)*(1.-a.f);

	        a.freq=a.MM+a.Mm/2.;
	        a.f=1.-a.Mm/(2*a.freq*(1-a.freq) );

#ifdef DEBUG
	       	std::cerr << "P:" << p << ", " << F << ", " << E << std::endl;
	       	std::cerr << iter <<" Pi= " <<  a.freq << ", Epsilon=" << a.error  << ", F=" << a.f << ", R=" << R[0]+R[1] << ": lnL=" << sumlnL << " : lnL="<< model.loglikelihood(site, a) << std::endl;
		std::cerr << ( pow(p,4)+pow(F,4)+pow(E,4) )/160000.0 << std::endl;
#endif
		a.f=F;
		a.freq=p;

		//BOUNDS CHECKS.
		
		float_t B= 5.1-float_t(iter)*0.015;

		if (fabs(R[0])>B){
			if(R[0]>0) R[0]=B;
			else R[0]=-B;
		}
		if (fabs(R[1])>B){
			if(R[1]>0) R[1]=B;
			else R[1]=-B;
		}

		a.freq-=R[0];
		a.f-=R[1];
	
        };

	a.error=0.75/(1.+exp(a.error) );
	a.freq=1./(1.+exp(a.freq) );
	a.f=1./(1.+exp(a.f) );

        a.MM=(1.-a.freq)*a.f;
        a.Mm=a.freq;
        a.mm=(1.-a.freq)*(1.-a.f);

	a.freq=a.MM+a.Mm/2.;
	a.f=1.-a.Mm/(2*a.freq*(1-a.freq) );

	a.coverage=site.getcoverage();

	count_t excluded=0;

        excluded=site.maskedcount();
	std::vector <float_t> temp_gofs(site.sample.size());
	a.gof=model.goodness_of_fit(site, a, temp_gofs, maxgof);
        a.efc=efc(site);

	if (iter==200) {
		std::cerr << "Failure to maximize " << iter << " " << a << "\n";
	}
	if ( a.gof<maxgof) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
                a.N-=1;
		return maximize_newton(site, a, model, gofs, maxgof, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return excluded;
}

count_t maximize_newton (Locus &site, Allele &a, models &model, std::vector <float_t> &gofs, const float_t &maxgof, const size_t &maxpitch)
{
	float_t J[3][3];		//The Jacobian/Hessian.
        float_t iJ[3][3]; 		//Inverse of the Jacobian/Hessian.
	float_t R[3]={100,100,100};	//The 'residual'. Th

        float_t det;

	std::vector <quartet_t>::const_iterator it=site.sample.begin();	//Lets us iterate over the quartets. 
	std::vector <quartet_t>::const_iterator end=site.sample.end();	 

	count_t iter=0;			//counts the number of iterations to let us know if we have a failure to converge.

	float_t MM=a.MM, Mm=a.Mm, mm=a.mm, p=a.freq, F=a.freq, E;//, q=1-a.freq;
//	a.f=log( (1.+p/(1.-p) )/(a.f+p/(1-p) )-1.);
//	a.f=log( (1.+(1.-p)/p )/(0.2+(1.-p)/p )-1.);
//	MM=pow(p, 2)+p*q*F;
//	Mm=2*p*q*(1-F);
//	mm=pow(q, 2)+p*q*F;

/*
	if (mm==0){
		if (Mm>MM) {mm+=0.01; Mm-=0.01;}
		else {mm+=0.01; MM-=0.01;}
	};
	if (Mm==0){
		if (mm>MM) {Mm+=0.01; mm-=0.01;}
		else {Mm+=0.01; MM-=0.01;}
	};
	if (MM==0){
		if (Mm>mm) {MM+=0.01; Mm-=0.01;}
		else {MM+=0.01; mm-=0.01;}
	};*/

	a.freq=Mm;
	a.f=MM/(MM+mm);


	if (a.freq==1.) a.freq-=0.000001;
	if (a.f==1.) a.f-=0.000001;
	if (a.freq==0.) a.freq+=0.000001;
	if (a.f==0.) a.f+=0.000001;
	if (a.error==0.5) a.error-=0.0000000001;
	if (a.error==0.0) a.error+=0.0000000001;

	a.f=log(1./a.f-1.);
	a.freq=log(1./a.freq-1.);
	a.error=log(0.75/a.error-1.);

	float_t deltalnL=10, maxlnL=DBL_MAX;

        while ( ( (fabs(R[0])+fabs(R[1])+fabs(R[2]) )>0.00001 || std::isnan(R[0]) || std::isnan(R[1]) || std::isnan(R[2]) ) && iter<200 && fabs(deltalnL)>0.0001 ){
#ifdef DEBUG
		std::cerr << fabs(R[0])+fabs(R[1])+fabs(R[2]) << ", " << iter << ", " << fabs(deltalnL) << std::endl;
#endif
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
		deltalnL=(maxlnL-sumlnL)+1;
		maxlnL=sumlnL;
/*		Check my matrix algebra.
		
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

		if (fabs(det)<0.00001) 
		{
			det < 0 ? det=-0.00001 : det=0.00001;
		}

		iJ[0][0]/=det; iJ[0][1]/=det; iJ[0][2]/=det;
		iJ[1][0]/=det; iJ[1][1]/=det; iJ[1][2]/=det;
		iJ[2][0]/=det; iJ[2][1]/=det; iJ[2][2]/=det;

                R[0]=(R[0]*iJ[0][0]+R[1]*iJ[0][1]+R[2]*iJ[0][2]);
                R[1]=(R[0]*iJ[1][0]+R[1]*iJ[1][1]+R[2]*iJ[1][2]);
                R[2]=(R[0]*iJ[2][0]+R[1]*iJ[2][1]+R[2]*iJ[2][2]);

		p=a.freq;
		F=a.f;
		E=a.error;

	        a.error=0.75/(1.+exp(a.error) );
	        a.freq=1./(1.+exp(a.freq) );
	        a.f=1./(1.+exp(a.f) );

	        a.MM=(1.-a.freq)*a.f;
	        a.Mm=a.freq;
	        a.mm=(1.-a.freq)*(1.-a.f);

	        a.freq=a.MM+a.Mm/2.;
	        a.f=1.-a.Mm/(2*a.freq*(1-a.freq) );

#ifdef DEBUG
	       	std::cerr << "P:" << p << ", " << F << ", " << E << std::endl;
	       	std::cerr << iter <<" Pi= " <<  a.freq << ", Epsilon=" << a.error  << ", F=" << a.f << ", R=" << R[0]+R[1]+R[2] << ": lnL=" << sumlnL << " : lnL="<< model.loglikelihood(site, a) << std::endl;
		std::cerr << ( pow(p,4)+pow(F,4)+pow(E,4) )/160000.0 << std::endl;
#endif
		a.f=F;
		a.freq=p;
		a.error=E;

		//BOUNDS CHECKS.
		
		float_t B= 5.1-float_t(iter)*0.015;

		if (fabs(R[0])>B){
			if(R[0]>0) R[0]=B;
			else R[0]=-B;
		}
		if (fabs(R[1])>B){
			if(R[1]>0) R[1]=B;
			else R[1]=-B;
		}
		if (fabs(R[2])>B){
			if(R[2]>0) R[2]=B;
			else R[2]=-B;
		}

		a.freq-=R[0];
		a.error-=R[1];
		a.f-=R[2];
	
        };

	a.error=0.75/(1.+exp(a.error) );
	a.freq=1./(1.+exp(a.freq) );
	a.f=1./(1.+exp(a.f) );

        a.MM=(1.-a.freq)*a.f;
        a.Mm=a.freq;
        a.mm=(1.-a.freq)*(1.-a.f);

	a.freq=a.MM+a.Mm/2.;
	a.f=1.-a.Mm/(2*a.freq*(1-a.freq) );

	a.coverage=site.getcoverage();

	count_t excluded=0;

        excluded=site.maskedcount();
	std::vector <float_t> temp_gofs(site.sample.size());
	a.gof=model.goodness_of_fit(site, a, temp_gofs, maxgof);
        a.efc=efc(site);

	if (iter==200) {
		std::cerr << "Failure to maximize " << iter << " " << a << "\n";
	}
	if ( a.gof<maxgof) {
		if (excluded==maxpitch){
			for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
			return excluded;
		}
                a.N-=1;
		return maximize_newton(site, a, model, gofs, maxgof, maxpitch);
	};
	for (size_t i=0; i<gofs.size(); i++) gofs[i]+=temp_gofs[i];
	return excluded;
}

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximize_grid (Locus &site, Allele &a, models &model, std::vector <float_t> &gofs, const float_t &MINGOF, const size_t &maxpitch)
{
	
	count_t N_=a.N;
	count_t P_=rint(a.MM*N_);
	count_t H_=rint(a.Mm*N_);
	count_t Q_=N_-P_-H_;
	if (P_+H_>N_) 
	{
		if(P_>N_)
		{
			P_=N_; 
			H_=0; 
			Q_=0;
		} else {
			Q_=0;
			H_=N_-P_;
		}
	}
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

/* Calculates the ratio of Major to minor allele, and reports chi-squared p-value of deviation from E=0.5.*/
void
get_bias (const Locus &site, Allele &a)
{
	std::vector <quartet_t>::const_iterator it=site.sample.begin(); 
	std::vector <quartet_t>::const_iterator end=site.sample.end(); 

	count_t tM=0, tm=0; 	//Total number of [M]ajor and [m]inor reads observed

	while (it!=end){
		if (!it->masked){
			if ( (*it)[a.major] != 0 && (*it)[a.minor] != 0) {tM+=(*it)[a.major]; tm+=(*it)[a.minor];}
		}
		++it;
	}

	a.bias=float_t(tM)/float_t(tM+tm);
	float_t half=float_t(tM+tm)/2.;
	float_t chi=powf(tM-half, 2)/half+powf(tm-half,2)/half;
	a.pbias=1.-gsl_cdf_chisq_P(chi, 1);
}

/* Uses a grid method to maximize the likelihood equations.*/
count_t maximize_analytical (Locus &site, Allele &a, models &model, std::vector <float_t> &gofs, const float_t &maxgof, const size_t &maxpitch)
{
	a.MM=1.0;
	a.Mm=0.0;
	a.mm=0.0;

	a.freq=1.0;

	//These are both undefined...
	a.f=0;
	a.gof=0;
	a.minor=4;
	float_t excluded=site.maskedcount();

	return excluded;
}
