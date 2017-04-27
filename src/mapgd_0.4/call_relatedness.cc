//relatedness from called sites.
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstring>
#include <sstream>

rel main (?){

	double mu, k1, k2, chisq;
	size_t N;
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (pow(N*2-1, 1), 2);
	gsl_matrix *Xm = gsl_matrix_alloc (pow(N*2-1, 1), 2);
	gsl_matrix *XmINV = gsl_matrix_alloc (pow(N*2-1, 1), 2);
	gsl_vector *yv = gsl_vector_alloc (pow(N*2-1, 1) );
	gsl_vector *wv = gsl_vector_alloc (pow(N*2-1, 1) );
	gsl_vector *cv = gsl_vector_alloc (2);
	gsl_matrix *cov = gsl_matrix_alloc (2, 2);
	{
		//relatedness to 0 ??
		for (size_t c2=1; c2<(N-1); c2++)
		{
				double d2=double(c2)/double(N*2-2);

				double mmmm=
				double Mmmm=
				double MMmm=
				double mmMm=
				double MmMm=
				double MMMm=
				double mmMM=
				double MmMM=
				double MMMM=

				double Den=mmmm+Mmmm+MMmm+mmMm+MmMm+MMMm+mmMM+MmMM+MMMM;

				double OX=((Mmmm+MmMm+MmMM)/2.+MMmm+MMMm+MMMM)/Den;
				double OY=((mmMm+MmMm+MMMm)/2.+mmMM+MmMM+MMMM)/Den;

				if (Den>10 && OX > 0 && OY > 0 && OX < 1 && OY <1 && OX!=0.5 && OY != 0.5)
				{
					double M=(OX+OY)/2.;

					Theta+=(mmmm*(0-OX)*(0-OY)+Mmmm*(1.0-OX)*(0-OY)/2.+MMmm*(1.-OX)*(0-OY)
								  +Mmmm*(0.0-OX)*(0-OY)/2.
						   +mmMm*(0-OX)*(1.0-OY)/2.+MmMm*(0.0-OX)*(0.0-OY)/4.+MMMm*(1.-OX)*(0.0-OY)/2.
						   +mmMm*(0-OX)*(0.0-OY)/2.+MmMm*(0.0-OX)*(1.0-OY)/4.+MMMm*(1.-OX)*(1.0-OY)/2.
						   			   +MmMm*(1.0-OX)*(0.0-OY)/4.
						   			   +MmMm*(1.0-OX)*(1.0-OY)/4.
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.0-OX)*(1.-OY)/2.+MMMM*(1.-OX)*(1.-OY) 
						   		       +MmMM*(1.0-OX)*(1.-OY)/2. )/sqrt(OX*(1-OX)*(1-OY)*OY );
					Theta_w+=Den;

					f_X+=(mmmm*(0-OX)*(0-OX)+Mmmm*(0.-OX)*(1.-OX)+MMmm*(1.-OX)*(1.-OX)
						   +mmMm*(0-OX)*(0-OX)+MmMm*(0.-OX)*(1.-OX)+MMMm*(1.-OX)*(1.-OX)
						   +mmMM*(0-OX)*(0-OX)+MmMM*(0.-OX)*(1.-OX)+MMMM*(1.-OX)*(1.-OX) )/(OX*(1-OX) );
					f_X_w+=Den;
					

					f_Y+=(mmmm*(0-OY)*(0-OY)+Mmmm*(0-OY)*(0-OY)+MMmm*(0-OY)*(0-OY)
						   +mmMm*(0.-OY)*(1.-OY)+MmMm*(0.-OY)*(1.-OY)+MMMm*(0.-OY)*(1.-OY)
						   +mmMM*(1.-OY)*(1.-OY)+MmMM*(1.-OY)*(1.-OY)+MMMM*(1.-OY)*(1.-OY) )/(OY*(1-OY ) );
					f_Y_w+=Den;
					if (M<0.45){
						gamma_XY+=(mmmm*(0-OX)*(0-OX)*(0-OY)+Mmmm*(0.-OX)*(1.-OX)*(0-OY)+MMmm*(1.-OX)*(1.-OX)*(0-OY)
								+(mmMm*(0-OX)*(0-OX)*(1.-OY)+MmMm*(1.-OX)*(0.-OX)*(1.-OY)+MMMm*(1.-OX)*(1.-OX)*(1.-OY) )/2.
								+(mmMm*(0-OX)*(0-OX)*(0.-OY)+MmMm*(1.-OX)*(0.-OX)*(0.-OY)+MMMm*(1.-OX)*(1.-OX)*(0.-OY) )/2.
								+mmMM*(0-OX)*(0-OX)*(1.-OY)+MmMM*(1.-OX)*(0.-OX)*(1.-OY)+MMMM*(1.-OX)*(1.-OX)*(1.-OY) )/(M*(1-M)*(1-2*M)); //p(1-p)(1-2p)
						gamma_XY_w+=Den;
	
						gamma_YX+=(mmmm*(0-OY)*(0-OY)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(0.0-OX)/2.+MMmm*(0.-OY)*(0.-OY)*(1.-OX)
										    +Mmmm*(0.-OY)*(0.-OY)*(1.0-OX)/2.
								+mmMm*(1.-OY)*(0.-OY)*(0.-OX)+MmMm*(1.-OY)*(0.-OY)*(0.0-OX)/2.+MMMm*(1.-OY)*(0.-OY)*(1.-OX)
										   	     +MmMm*(1.-OY)*(0.-OY)*(1.0-OX)/2.
								+mmMM*(1.-OY)*(1.-OY)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(0.0-OX)/2.+MMMM*(1.-OY)*(1.-OY)*(1.-OX) 
										    	     +MmMM*(1.-OY)*(1.-OY)*(1.0-OX)/2.
														)/(M*(1-M)*(1-2*M) );
						gamma_YX_w+=Den;
					}
					mu=(mmmm*(0-OY)*(0-OY)*(0-OX)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(1.-OX)*(0.0-OX)+MMmm*(0.-OY)*(0.-OY)*(1.-OX)*(1.-OX)
							+mmMm*(1.-OY)*(0.-OY)*(0.-OX)*(0.-OX)+MmMm*(0.-OY)*(1.-OY)*(0.-OX)*(1.-OX)+MMMm*(0.-OY)*(1.-OY)*(1.-OX)*(1.-OX)
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(1.-OX)*(0.-OX)+MMMM*(1.-OY)*(1.-OY)*(1.-OX)*(1.-OX) )/Den;
					double k1x=pow(OX*(1-OX), 2);
					double k2x=OX*(1-OX)*(3*OX*OX-3*OX+1);
					double k1y=pow(OY*(1-OY), 2);
					double k2y=OY*(1-OY)*(3*OY*OY-3*OY+1);
					k1=sqrt(k1x*k1y);
					k2=sqrt(k2x*k2y);
					gsl_matrix_set (Xm, c2, 0, k1);
					gsl_matrix_set (Xm, c2, 1, k2);
      
					gsl_vector_set (yv, c2, mu);
					gsl_vector_set (wv, c2, Den);
			}
		}
		#define beta(i) (gsl_vector_get(cv,(i)))
		gsl_multifit_wlinear (Xm, wv, yv, cv, cov, &chisq, work);
		std::cout << X << "\t" << Y << "\t" << f_X/f_X_w << "\t" << f_X_w << "\t" << f_Y/f_Y_w << "\t" << f_Y_w << "\t" << Theta/Theta_w << "\t" <<Theta_w <<"\t" << gamma_XY/gamma_XY_w << "\t" << gamma_XY_w << "\t" << gamma_YX/gamma_YX_w << "\t" << gamma_YX_w << "\t" << beta(1) << "\t" << "0" << "\t" << beta(0) << "\t" << std::endl;

	}
}	
}

/* rho_xy|z=rho_xy-rho_xz*rho_yz/( sqrt(1-rho_xz^2)*rho(1-rho_yz^2) */
	std::cout << "@END_TABE\n";
}	
