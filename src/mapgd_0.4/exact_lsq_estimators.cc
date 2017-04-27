#include "circular_list.h"
#include <gsl/gsl_multifit.h>

/*These are the beginings of a least squared method. None of this is used for anything, just to remind me of the equations*/

class Slow_index
{
	private:
	size_t _1, _2, _3, _4, _5;
	public:
	Slow_index(const size_t &_N)
	{	
		_1=9*_N*_N*_N; 
		_2=9*_N*_N;
		_3=9*_N; 
		_4=9;
		_5=3;
	}

	//d1(N), x(N), y(N), i(3), j(3) 

	const size_t 	
	operator () (const uint32_t &c1, const uint32_t &c2, const uint32_t &x, const uint32_t &y, const uint32_t &i, const uint32_t &j)	{return c1*_1+c2*_2+x*_3+y*_4+i*_5+j;};
};

static uint32_t Mask[32]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
			0x00000010, 0x00000020, 0x00000040, 0x00000080,
			0x00000100, 0x00000200, 0x00000400, 0x00000800,
			0x00001000, 0x00002000, 0x00004000, 0x00008000,
			0x00010000, 0x00020000, 0x00040000, 0x00080000,
			0x00100000, 0x00200000, 0x00400000, 0x00800000,
			0x01000000, 0x02000000, 0x04000000, 0x08000000,
			0x10000000, 0x20000000, 0x40000000, 0x80000000};

inline void fast_bit_sumN(const uint32_t *place, const uint32_t &N, const uint32_t &bit, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
	uint32_t mask=Mask[bit];
	while (this_place!=end) {
		*sum+=*this_place & mask;
		this_place++;
	}
}

inline int get(const uint32_t *place, const uint32_t &x, const uint32_t &bit)
{
	//nstd::cerr << "(" << x << ")" << std::endl;
	if ((*(place+x) & Mask[bit])!=0 ){
		return 1;
	} else  {
		return 0;
	}
}

inline void fast_het_sumN(const uint32_t *place, const uint32_t &N, const uint32_t &bit, uint32_t *sum)
{
	const uint32_t *end=place+N;
	const uint32_t *this_place=place;
	uint32_t mask=Mask[bit];
	while (this_place!=end) {
		//if((*this_place & mask)==0) std::cerr << 0; 
		//else std::cerr << 1;
		//if((*(this_place+1) & mask)==0) std::cerr << 0; 
		//else std::cerr << 1;
		*sum+=( (*this_place & mask)!=(*(++this_place) & mask) );
		this_place++;
	}
}

/*const size_t & fast_index(const uint &d, const uint32_t &i, const uint32_t &j)
{
	return _index[d][i][j];
}*/

int main (int argc, char **argv){

uint32_t N=atoi(argv[1]);

uint32_t SIZE=9*N*N*N*N; 
uint32_t BLOCK_SIZE=N;
uint32_t *D=new uint32_t [SIZE];

uint32_t * bit=new uint32_t [N];

Slow_index slow_index(N);

std::istream *in=&std::cin;
	
uint32_t readed=0;

while(in->read((char *)(bit), BLOCK_SIZE*sizeof(uint32_t) ) )
{
	for (size_t x=0; x<N; x+=2)
	{
		for (size_t y=x+1; y<N; y+=2)
		{
			for (size_t k=0; k<32; k++)
			{
				int i=(get(bit, x, k)+get(bit, x+1, k) );
				int j=(get(bit, y, k)+get(bit, y+1, k) );
				uint32_t d1=0;
				uint32_t d2=0;
				fast_bit_sumN(bit, N, k, &d1);
				d2=d1;
				d1-=i;
				d2-=j;
				(*(D+slow_index(d1, d2, x, y, i, j) ) )++;
			}
		}
	}
}

std::cout << "x y f_X f_Y Theta gamma_XY gamma_YX Delta delta\n";

/* rho_xy|z=rho_xy-rho_xz*rho_yz/( sqrt(1-rho_xz^2)*rho(1-rho_yz^2) */

double f_X=0, f_X_w=0, f_Y=0, f_Y_w=0, Theta=0, Theta_w=0, gamma_XY=0, gamma_XY_w=0, gamma_YX=0, gamma_YX_w=0, Z=0, Z_w=0;
double mu, k1, k2;


//    gsl_multifit_linear_workspace * work 
//      = gsl_multifit_linear_alloc (n, 3);

gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (pow(N-1, 2 ), 2);
gsl_matrix *X = gsl_matrix_alloc (pow(N-1,2), 2);
gsl_vector *y = gsl_vector_alloc (pow(N-1,2) );
gsl_vector *w = gsl_vector_alloc (pow(N-1,2) );
gsl_vector *c = gsl_vector_alloc (3);
gsl_matrix *cov = gsl_matrix_alloc (3, 3);
double chisq;


for (size_t x=0; x<N/2; x++)
{
	for (size_t Y=x+1; Y<N/2; Y++)
	{
		for (size_t c1=1; c1<N-1; c1++)
		{
			for (size_t c2=1; c2<N-1; c2++)
			{
				double d1=double(c1)/double(N/2-2);
				double d2=double(c2)/double(N/2-2);

				double mmmm=double(*(D+slow_index(c1, c2, x, Y, 0, 0)));
				double Mmmm=double(*(D+slow_index(c1, c2, x, Y, 1, 0)));
				double MMmm=double(*(D+slow_index(c1, c2, x, Y, 2, 0)));
				double mmMm=double(*(D+slow_index(c1, c2, x, Y, 0, 1)));
				double MmMm=double(*(D+slow_index(c1, c2, x, Y, 1, 1)));
				double MMMm=double(*(D+slow_index(c1, c2, x, Y, 2, 1)));
				double mmMM=double(*(D+slow_index(c1, c2, x, Y, 0, 2)));
				double MmMM=double(*(D+slow_index(c1, c2, x, Y, 1, 2)));
				double MMMM=double(*(D+slow_index(c1, c2, x, Y, 2, 2)));

				double Den=mmmm+Mmmm+MMmm+mmMm+MmMm+MMMm+mmMM+MmMM+MMMM;

				if (Den>10)
				{
					double OX=(Mmmm+MmMm+MmMM)/2.+MMmm+MMMm+MMMM;
					double OY=(mmMm+MmMm+MMMm)/2.+mmMM+MmMM+MMMM;
					double M=(OX+OY)/2.;

					Theta+=Den*(mmmm*double(0-OX)*double(0-OY)+Mmmm*double(0.5-OX)*double(0-OY)+MMmm*double(1.-OX)*double(0-OY)
						   +mmMm*(0-OX)*(0.5-OY)+MmMm*(0.5-OX)*(0.5-OY)+MMMm*(1.-OX)*(0.5-OY)
						   +mmMM*(0-OX)*(1.-OY)+MmMM*(0.5-OX)*(1.-OY)+MMMM*(1.-OX)*(1.-OY) )/sqrt(OX*(1-OX)*(1-OY)*OY );
					Theta_w+=Den;

					f_X+=Den*(mmmm*double(0-OX)*double(0-OX)+Mmmm*double(0.5-OX)*double(0.5-OX)+MMmm*double(1.-OX)*double(1.-OX)
						   +mmMm*(0-OX)*(0-OX)+MmMm*(0.5-OX)*(0.5-OX)+MMMm*(1.-OX)*(1.-OX)
						   +mmMM*(0-OX)*(0-OX)+MmMM*(0.5-OX)*(0.5-OX)+MMMM*(1.-OX)*(1.-OX) )/(OX*(1-OX) );
					f_X_w+=Den;
					

					f_Y+=Den*(mmmm*double(0-OY)*double(0-OY)+Mmmm*double(0-OY)*double(0-OY)+MMmm*double(0-OY)*double(0-OY)
						   +mmMm*double(0.5-OY)*double(0.5-OY)+MmMm*double(0.5-OY)*double(0.5-OY)+MMMm*double(0.5-OY)*double(0.5-OY)
						   +mmMM*(1.-OY)*(1.-OY)+MmMM*(1.-OY)*(1.-OY)+MMMM*(1.-OY)*(1.-OY) )/(OY*(1-OY ) );
					f_Y_w+=Den;

					gamma_XY+=Den*(mmmm*(0-OX)*(0-OX)*(0-OY)+Mmmm*(0.5-OX)*(0.5-OX)*(0-OY)+MMmm*(1.-OX)*(1.-OX)*(0-OY)
							+mmMm*(0-OX)*(0-OX)*(0.5-OY)+MmMm*(0.5-OX)*(0.5-OX)*(0.5-OY)+MMMm*(1.-OX)*(1.-OX)*(0.5-OY)
							+mmMM*(0-OX)*(0-OX)*(1.-OY)+MmMM*(0.5-OX)*(0.5-OX)*(1.-OY)+MMMM*(1.-OX)*(1.-OX)*(1.-OY) )/sqrt( (1-2*OX)/sqrt(OX*(1-OX) )*(1-2*OY)/sqrt(OY*(1-OY) ) );
					gamma_XY_w+=Den;

					gamma_YX+=Den*(mmmm*(0-OY)*(0-OY)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(0.5-OX)+MMmm*(0.-OY)*(0.-OY)*(1.-OX)
							+mmMm*(0.5-OY)*(0.5-OY)*(0.-OX)+MmMm*(0.5-OY)*(0.5-OY)*(0.5-OX)+MMMm*(0.5-OY)*(0.5-OY)*(1.-OX)
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(0.5-OX)+MMMM*(1.-OY)*(1.-OY)*(1.-OX) )/(sqrt( (1-2*OX)/sqrt(OX*(1-OX) )*(1-2*OY)/sqrt(OY*(1-OY) ) ) );
					gamma_YX_w+=Den;

					mu=(mmmm*(0-OY)*(0-OY)*(0-OX)*(0-OX)+Mmmm*(0.-OY)*(0.-OY)*(0.5-OX)*(0.5-OX)+MMmm*(0.-OY)*(0.-OY)*(1.-OX)*(1.-OX)
							+mmMm*(0.5-OY)*(0.5-OY)*(0.-OX)*(0.-OX)+MmMm*(0.5-OY)*(0.5-OY)*(0.5-OX)*(0.5-OX)+MMMm*(0.5-OY)*(0.5-OY)*(1.-OX)*(1.-OX)
							+mmMM*(1.-OY)*(1.-OY)*(0.-OX)*(0.-OX)+MmMM*(1.-OY)*(1.-OY)*(0.5-OX)*(0.5-OX)+MMMM*(1.-OY)*(1.-OY)*(1.-OX)*(1.-OX) );
					k1=pow(M*(1-M) ,2);
					k2=1/(1-M)+1/M-3;

					gsl_matrix_set (X, c1*(N-1)+c2, 0, k1);
					gsl_matrix_set (X, c1*(N-1)+c2, 1, k2);
      
					gsl_vector_set (y, c1*(N-1)+c2, mu);
					gsl_vector_set (w, c1*(N-1)+c2, Den);
				}
			}
		}
		gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);
		gsl_multifit_linear_free (work);

#define beta(i) (gsl_vector_get(c,(i)))

		std::cout << x << " " << Y << " " << f_X/f_X_w << " " << f_Y/f_Y_w << " " << Theta/Theta_w <<" " << gamma_XY/gamma_XY_w << " " << gamma_YX/gamma_YX_w << " " << beta(1) << " " << beta(2) << std::endl;
	}
}	
}	
