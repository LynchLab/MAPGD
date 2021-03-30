#include "relatedness.h"

#define BUFFER_SIZE 500

#ifndef NOGSL

///Takes a thing and does something.
size_t 
freqtoi(float_t in)
{
	return size_t(in*E_LIM*2) < E_LIM ? size_t(in*E_LIM*2) : E_LIM-1;
}

#define ITER_MAX    150

double
rel_ll2 (const Relatedness &rel, std::vector < std::pair<Genotype_pair_tuple, size_t> > *hashed_genotypes_p)
{
    gsl_vector *v=gsl_vector_alloc(7);

	gsl_vector_set(v, 0, rel.f_X_);
	gsl_vector_set(v, 1, rel.f_Y_);
	gsl_vector_set(v, 2, rel.theta_XY_);
	gsl_vector_set(v, 3, rel.gamma_XY_);
	gsl_vector_set(v, 4, rel.gamma_YX_);
	gsl_vector_set(v, 5, rel.Delta_XY_);
	gsl_vector_set(v, 6, rel.delta_XY_);
	
	double ret=rel_ll(v, (void *)(hashed_genotypes_p) );
    gsl_vector_free(v);
    return ret;
}


int
get_invJ(gsl_matrix *invJ, const Relatedness &a, std::map <Genotype_pair_tuple, size_t> &counts)
{
	gsl_matrix *J=gsl_matrix_alloc(7,7);

    float_t (*Jptr[7][7]) (const Genotype_pair &, const Constants <float_t, const std::pair <const Genotype_pair &, const Relatedness &> > &);
        
    Jptr[0][0]=&J00; Jptr[0][1]=&J01; Jptr[0][2]=&J02; Jptr[0][3]=&J03; Jptr[0][4]=&J04; Jptr[0][5]=&J05; Jptr[0][6]=&J06;
    Jptr[1][0]=&J10; Jptr[1][1]=&J11; Jptr[1][2]=&J12; Jptr[1][3]=&J13; Jptr[1][4]=&J14; Jptr[1][5]=&J15; Jptr[1][6]=&J16;
    Jptr[2][0]=&J20; Jptr[2][1]=&J21; Jptr[2][2]=&J22; Jptr[2][3]=&J23; Jptr[2][4]=&J24; Jptr[2][5]=&J25; Jptr[2][6]=&J26;
    Jptr[3][0]=&J30; Jptr[3][1]=&J31; Jptr[3][2]=&J32; Jptr[3][3]=&J33; Jptr[3][4]=&J34; Jptr[3][5]=&J35; Jptr[3][6]=&J36;
    Jptr[4][0]=&J40; Jptr[4][1]=&J41; Jptr[4][2]=&J42; Jptr[4][3]=&J43; Jptr[4][4]=&J44; Jptr[4][5]=&J45; Jptr[4][6]=&J46;
    Jptr[5][0]=&J50; Jptr[5][1]=&J51; Jptr[5][2]=&J52; Jptr[5][3]=&J53; Jptr[5][4]=&J54; Jptr[5][5]=&J55; Jptr[5][6]=&J56;
    Jptr[6][0]=&J60; Jptr[6][1]=&J61; Jptr[6][2]=&J62; Jptr[6][3]=&J63; Jptr[6][4]=&J64; Jptr[6][5]=&J65; Jptr[6][6]=&J66;

    std::map <Genotype_pair_tuple, size_t>::iterator it, end;
	Genotype_pair v;
	Constants <float_t, const std::pair<const Genotype_pair &, const Relatedness &> > consts(REL_CNTS, REL_ARRAY);
    std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(counts.begin(), counts.end() );
    
	size_t c;

	float_t last_m=0;

	it=counts.begin();
	end=counts.end();

    gsl_matrix_set_zero(J);

	while (it!=end ){
		v=Genotype_pair::from_tuple(it->first);
		c=it->second;
		if (v.m>0.0)
		{
			if (v.m!=last_m){
				consts.recalculate( std::pair <Genotype_pair, Relatedness> (v, a) );
				last_m=v.m;
			}
			for (int x=0; x<7; x++) {
				for (int y=0; y<7; y++) {
					gsl_matrix_set(J, x, y, gsl_matrix_get(J,x,y)-(*Jptr[x][y])(v, consts)*c);
				}
			}
        }
        it++;
    }
    int signum=0;
	gsl_permutation * p = gsl_permutation_alloc (7);
    gsl_linalg_LU_decomp(J,p, &signum);
    int ret=gsl_linalg_LU_invert(J, p, invJ);
	gsl_permutation_free (p);
    gsl_matrix_free (J);
    return ret;
}

int 
newton (Relatedness &a, std::map <Genotype_pair_tuple, size_t> &counts)
{

    float_t (*Jptr[7][7]) (const Genotype_pair &, const Constants <float_t, const std::pair <const Genotype_pair &, const Relatedness &> > &);
        
    Jptr[0][0]=&J00; Jptr[0][1]=&J01; Jptr[0][2]=&J02; Jptr[0][3]=&J03; Jptr[0][4]=&J04; Jptr[0][5]=&J05; Jptr[0][6]=&J06;
    Jptr[1][0]=&J10; Jptr[1][1]=&J11; Jptr[1][2]=&J12; Jptr[1][3]=&J13; Jptr[1][4]=&J14; Jptr[1][5]=&J15; Jptr[1][6]=&J16;
    Jptr[2][0]=&J20; Jptr[2][1]=&J21; Jptr[2][2]=&J22; Jptr[2][3]=&J23; Jptr[2][4]=&J24; Jptr[2][5]=&J25; Jptr[2][6]=&J26;
    Jptr[3][0]=&J30; Jptr[3][1]=&J31; Jptr[3][2]=&J32; Jptr[3][3]=&J33; Jptr[3][4]=&J34; Jptr[3][5]=&J35; Jptr[3][6]=&J36;
    Jptr[4][0]=&J40; Jptr[4][1]=&J41; Jptr[4][2]=&J42; Jptr[4][3]=&J43; Jptr[4][4]=&J44; Jptr[4][5]=&J45; Jptr[4][6]=&J46;
    Jptr[5][0]=&J50; Jptr[5][1]=&J51; Jptr[5][2]=&J52; Jptr[5][3]=&J53; Jptr[5][4]=&J54; Jptr[5][5]=&J55; Jptr[5][6]=&J56;
    Jptr[6][0]=&J60; Jptr[6][1]=&J61; Jptr[6][2]=&J62; Jptr[6][3]=&J63; Jptr[6][4]=&J64; Jptr[6][5]=&J65; Jptr[6][6]=&J66;

	gsl_matrix *J=gsl_matrix_alloc(7,7);
	gsl_vector *R=gsl_vector_alloc(7);
	gsl_vector *new_R=gsl_vector_alloc(7);
	gsl_permutation * p = gsl_permutation_alloc (7);

	std::map <Genotype_pair_tuple, size_t>::iterator it, end;
	Genotype_pair v;
	Constants <float_t, const std::pair<const Genotype_pair &, const Relatedness &> > consts(REL_CNTS, REL_ARRAY);
    std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(counts.begin(), counts.end() );

    /*
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_real_distribution<double> unif(-0.5, 1.);
    */

  //  std::cerr << "HI!\n";

	size_t c;

	float_t last_m=0;
	int iter=0;
	int size=counts.size();
    int I[7]={0};	
	while (true){
//    std::cerr << "HI!2\n";
		it=counts.begin();
		end=counts.end();
		gsl_matrix_set_zero(J);
		gsl_vector_set_zero(R);

	while (it!=end ){
		v=Genotype_pair::from_tuple(it->first);
		c=it->second;
		if (v.m>0.0)
		{
			if (v.m!=last_m){
				consts.recalculate( std::pair <Genotype_pair, Relatedness> (v, a) );
				last_m=v.m;
			}
			for (int x=0; x<7; x++) {
				for (int y=0; y<7; y++) {
					gsl_matrix_set(J, x, y, gsl_matrix_get(J,x,y)+(*Jptr[x][y])(v, consts)*c);
				}
			}

			gsl_vector_set(R, 0, gsl_vector_get(R, 0)-H0(v, consts)*c );
			gsl_vector_set(R, 1, gsl_vector_get(R, 1)-H1(v, consts)*c );
			gsl_vector_set(R, 2, gsl_vector_get(R, 2)-H2(v, consts)*c );
			gsl_vector_set(R, 3, gsl_vector_get(R, 3)-H3(v, consts)*c );
			gsl_vector_set(R, 4, gsl_vector_get(R, 4)-H4(v, consts)*c );
			gsl_vector_set(R, 5, gsl_vector_get(R, 5)-H5(v, consts)*c );
			gsl_vector_set(R, 6, gsl_vector_get(R, 6)-H6(v, consts)*c );

           //std::cerr << 
	
		}
		++it;
	//#pragma omp taskwait
	}
	if (gsl_blas_dnrm2(R)< 0.0001) {
    //    std::cerr << "success!!\n";
        break;
	}
	if (std::isnan(gsl_blas_dnrm2(R))){
    //    std::cerr << "NAN!!\n";
        break;
	}
	if (gsl_matrix_isnull(J)){
    //    std::cerr << "NULL!!\n";
        break;
	}
    //printf("iter? %d\n", iter);
	if (iter>=ITER_MAX){
     //   std::cerr << "Itermax reached\n";
        break;
	}
	++iter;

	int signum=0;
    
    gsl_linalg_LU_decomp (J, p, &signum);
    gsl_set_error_handler_off();
    /*printf( "J:%d, %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", iter, gsl_matrix_get(J, 0, 0), gsl_matrix_get(J,1,0), gsl_matrix_get(J,2,0), gsl_matrix_get(J,3,0), gsl_matrix_get(J,4,0), gsl_matrix_get(J,5,0), gsl_matrix_get(J,6,0) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 1), gsl_matrix_get(J,1,1), gsl_matrix_get(J,2,1), gsl_matrix_get(J,3,1), gsl_matrix_get(J,4,1), gsl_matrix_get(J,5,1), gsl_matrix_get(J,6,1) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 2), gsl_matrix_get(J,1,2), gsl_matrix_get(J,2,2), gsl_matrix_get(J,3,2), gsl_matrix_get(J,4,2), gsl_matrix_get(J,5,2), gsl_matrix_get(J,6,2) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 3), gsl_matrix_get(J,1,3), gsl_matrix_get(J,2,3), gsl_matrix_get(J,3,3), gsl_matrix_get(J,4,3), gsl_matrix_get(J,5,3), gsl_matrix_get(J,6,3) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 4), gsl_matrix_get(J,1,4), gsl_matrix_get(J,2,4), gsl_matrix_get(J,3,4), gsl_matrix_get(J,4,4), gsl_matrix_get(J,5,4), gsl_matrix_get(J,6,4) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 5), gsl_matrix_get(J,1,5), gsl_matrix_get(J,2,5), gsl_matrix_get(J,3,5), gsl_matrix_get(J,4,5), gsl_matrix_get(J,5,5), gsl_matrix_get(J,6,5) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 6), gsl_matrix_get(J,1,6), gsl_matrix_get(J,2,6), gsl_matrix_get(J,3,6), gsl_matrix_get(J,4,6), gsl_matrix_get(J,5,6), gsl_matrix_get(J,6,6) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 7), gsl_matrix_get(J,1,7), gsl_matrix_get(J,2,7), gsl_matrix_get(J,3,7), gsl_matrix_get(J,4,7), gsl_matrix_get(J,5,7), gsl_matrix_get(J,6,7) );
    printf( "      %.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n", gsl_matrix_get(J, 0, 8), gsl_matrix_get(J,1,8), gsl_matrix_get(J,2,8), gsl_matrix_get(J,3,8), gsl_matrix_get(J,4,8), gsl_matrix_get(J,5,8), gsl_matrix_get(J,6,8) );
	*/
    int status=gsl_linalg_LU_solve (J, p, R, new_R);
    if (status!=0) {
	    gsl_permutation_free (p);
    	gsl_vector_free (R);
    	gsl_vector_free (new_R);
    	gsl_matrix_free (J);
        std::cerr << "WTF?!?\n";
        return 1;
    }
    gsl_vector_memcpy (R, new_R);

    double k[7];
    for(int i=0; i<7; i++) k[i]=1.-1./(1.+exp(float(++I[i])/2.-2.) );

    //printf( "T:%d, %.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", iter, a.f_X_, a.f_Y_, a.theta_XY_, a.gamma_XY_, a.gamma_YX_, a.Delta_XY_, a.delta_XY_);
    //printf( "R:%d, %.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", iter, gsl_vector_get(R, 0), gsl_vector_get(R,1), gsl_vector_get(R,2), gsl_vector_get(R,3), gsl_vector_get(R,4), gsl_vector_get(R,5), gsl_vector_get(R,6) );

    a.f_X_+=gsl_vector_get(R, 0)*k[0];
	a.f_Y_+=gsl_vector_get(R,1)*k[1];
	a.theta_XY_+=gsl_vector_get(R,2)*k[2];
	a.gamma_XY_+=gsl_vector_get(R,3)*k[3];
	a.gamma_YX_+=gsl_vector_get(R,4)*k[4];
	a.delta_XY_+=gsl_vector_get(R,5)*k[5];
	a.Delta_XY_+=gsl_vector_get(R,6)*k[6];

    while(std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) )
    {
        //printf( "I:%d, %.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", iter, a.f_X_, a.f_Y_, a.theta_XY_, a.gamma_XY_, a.gamma_YX_, a.Delta_XY_, a.delta_XY_);

        //std::cerr << iter << ", " << k << std::endl;
	    a.f_X_-=gsl_vector_get(R, 0)*k[0];
        a.f_Y_-=gsl_vector_get(R,1)*k[1];
        a.theta_XY_-=gsl_vector_get(R,2)*k[2];
        a.gamma_XY_-=gsl_vector_get(R,3)*k[3];
        a.gamma_YX_-=gsl_vector_get(R,4)*k[4];
        a.delta_XY_-=gsl_vector_get(R,5)*k[5];
        a.Delta_XY_-=gsl_vector_get(R,6)*k[6];
    
        for(int i=0; i<7; i++) k[i]=1.-1./(1.+exp(float(I[i]-=2)/2.-2.) );
    
    /*    int i=0;
        a.f_X_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.f_X_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=1;
        a.f_Y_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.f_Y_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=2;
        a.theta_XY_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.theta_XY_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=3;
        a.gamma_XY_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.gamma_XY_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=4;
        a.gamma_YX_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.gamma_YX_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=5;
        a.Delta_XY_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.Delta_XY_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );

        i=6;
        a.delta_XY_+=gsl_vector_get(R, i)*k[i];
        if( std::isnan(rel_ll2(a, &hashed_genotypes_vector) ) ) I[i]-=3;
        a.delta_XY_-=gsl_vector_get(R, i)*k[i];
        k[i]=1.-1./(1.+exp(float(I[i])/2.-2.) );*/

        for (int j=0; j<7; j++){
          // printf("%d,", I[j]);
           if (I[j]<-30) {
               iter=ITER_MAX;
               return 1;
           }
        }
        //printf("\n");

	    a.f_X_+=gsl_vector_get(R, 0)*k[0];//+uinf(eng);
        a.f_Y_+=gsl_vector_get(R, 1)*k[1];//+uinf(eng);
        a.theta_XY_+=gsl_vector_get(R,2)*k[2];//+unif(eng);
        a.gamma_XY_+=gsl_vector_get(R,3)*k[3];//+unif(eng);
        a.gamma_YX_+=gsl_vector_get(R,4)*k[4];//+unif(eng);
        a.delta_XY_+=gsl_vector_get(R,5)*k[5];//+unif(eng);
        a.Delta_XY_+=gsl_vector_get(R,6)*k[6];//+unif(eng);

/*        memset(&I, 0, sizeof(int)*7);
	    a.f_X_=unif(eng);
        a.f_Y_=unif(eng);
        a.theta_XY_=unif(eng);
        a.gamma_XY_=unif(eng);
        a.gamma_YX_=unif(eng);
        a.delta_XY_=unif(eng);
        a.Delta_XY_=unif(eng);
        //printf( "%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", a.f_X_, a.f_Y_, a.theta_XY_, a.gamma_XY_, a.gamma_YX_, a.Delta_XY_, a.delta_XY_);
 */   }
    //printf("out!\n");


	}
	gsl_permutation_free (p);
	gsl_vector_free (R);
	gsl_vector_free (new_R);
	gsl_matrix_free (J);

    if ( iter>=ITER_MAX) return 1;
	if ( fabs(a.f_X_)>1 or fabs(a.f_Y_)>1  or fabs(a.theta_XY_) >1 or fabs(a.gamma_XY_)>1 or fabs(a.gamma_YX_)>1 or fabs(a.delta_XY_)>1 or fabs(a.Delta_XY_)>1 ) return 1; 
    return 0;
}

std::map <Genotype_pair_tuple, size_t> 
hash_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y, const bool &l2o, const bool &called)
{
	std::stringstream fb_copy(file_buffer.str() );

    Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	Population genotypes=gcf_in.read_header();
	float_t N=genotypes.likelihoods.size();
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		float s=(genotypes.likelihoods[x].mm+genotypes.likelihoods[y].mm)*2.+(genotypes.likelihoods[x].Mm+genotypes.likelihoods[y].Mm);
		if (genotypes.likelihoods[x].N>1 && genotypes.likelihoods[y].N>1 && genotypes.m>0 && genotypes.m<1  ){
            if (called) {
                if (l2o) counts[convert_called(genotypes.likelihoods[x], genotypes.likelihoods[y], (2.*genotypes.m*N-s)/(2.*N-4.), 3)]+=1;
    			else counts[convert_called(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 3)]+=1;
            } else {
                if (l2o) counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], (2.*genotypes.m*N-s)/(2.*N-4.), 3)]+=1;
    			else counts[convert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 3)]+=1;
            }
        }
	}
	return counts;
}


//TODO This is a bottle neck, surly we can make if faster.
std::map <Genotype_pair_tuple, size_t> 
downsample_genotypes (const std::stringstream &file_buffer, const size_t &x, const size_t &y, const bool &l2o, const bool &called)
{
	std::stringstream fb_copy(file_buffer.str() );
	
	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	gcf_in.open(&fb_copy, std::ios::in );
	Population genotypes=gcf_in.read_header();
	float_t N=genotypes.likelihoods.size();
	std::map <Genotype_pair_tuple, size_t> counts;
	while(gcf_in.table_is_open() ){
		gcf_in.read(genotypes);
		float s=(genotypes.likelihoods[x].mm+genotypes.likelihoods[y].mm)*2.+(genotypes.likelihoods[x].Mm+genotypes.likelihoods[y].Mm);
		if (genotypes.likelihoods[x].N>1 && genotypes.likelihoods[y].N>1 && genotypes.m>0 && genotypes.m<1  ){
			if (l2o) counts[downvert(genotypes.likelihoods[x], genotypes.likelihoods[y], (2.*genotypes.m*N-s)/(2.*N-4.), 3)]+=1;
			else counts[downvert(genotypes.likelihoods[x], genotypes.likelihoods[y], genotypes.m, 3)]+=1;
		}
	}
	return counts;
}

//Relatedness global_relatedness;

/*Does a regression of allele frequency of the samples on the population allele frequency*/
void 
set_e(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes)
{
	std::map<Genotype_pair_tuple, size_t>::iterator start=hashed_genotypes.begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=hashed_genotypes.end();
	std::map<Genotype_pair_tuple, size_t>::iterator it=start;

	Genotype_pair pair;

/*
        float_t X_MM;
        float_t X_Mm;
        float_t X_mm;
        float_t Y_MM;
        float_t Y_Mm;
        float_t Y_mm;
        float_t m;
*/
	float_t freq[E_LIM]={0}, sum[E_LIM]={0};
        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
		freq[freqtoi(pair.m)]+=pair.m;
		sum[freqtoi(pair.m)]+=it->second;
		relatedness.e_X_[freqtoi(pair.m)]+=exp(-pair.X_mm)+exp(-pair.X_Mm)/2.;
		relatedness.e_Y_[freqtoi(pair.m)]+=exp(-pair.Y_mm)+exp(-pair.Y_Mm)/2.;
		++it;
	}

	for (size_t x=0; x<E_LIM; ++x){
		relatedness.e_X_[x]=(relatedness.e_X_[x]-freq[x])/sum[x];
		relatedness.e_Y_[x]=(relatedness.e_Y_[x]-freq[x])/sum[x];
//		std::cerr << x << ": " << relatedness.e_X_[x] << ", " << relatedness.e_Y_[x] << std::endl;
	}
	//global_relatedness=relatedness;
}


/*Guess starting values of relatedness for the maximization procedure*/
void 
gestimate(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &counts)
{
	relatedness.zero();
    relatedness.f_X_=0.3333;
    relatedness.f_Y_=0.3333;
    relatedness.gamma_XY_=0.8333;
    relatedness.gamma_YX_=0.8333;
    relatedness.delta_XY_=0.1111;
    /*	std::map<Genotype_pair_tuple, size_t>::iterator start=counts.begin();
	std::map<Genotype_pair_tuple, size_t>::iterator end=counts.end();
	std::map<Genotype_pair_tuple, size_t>::iterator it=start;
	Genotype_pair pair;

	for (size_t x=0; x<?; x++){
	?
	double OX=
	double OY=

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
#define beta(i) (gsl_vector_get(cv,(i)))
	gsl_multifit_wlinear (Xm, wv, yv, cv, cov, &chisq, work);
        std::cerr << "\t" << f_X_w << "\t" << f_Y/f_Y_w << "\t" << f_Y_w << "\t" << Theta/Theta_w << "\t" <<Theta_w <<"\t" << gamma_XY/gamma_XY_w << "\t" << gamma_XY_w << "\t" << gamma_YX/gamma_YX_w << "\t" << gamma_YX_w << "\t" << beta(1) << "\t" << "0" << "\t" << beta(0) << "\t" << std::endl;
*/
/*        while(it!=end){
		pair=Genotype_pair::from_tuple(it->first);
               	inc_guess(relatedness, pair, it->second);
                ++(it);
                ++(it);
        }
*/

}

void 
inc_f(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	float_t P2=P*P;
	float_t var=P-P2;
	float_t denom=pow(var, 2);
	if(pair.m!=0){
		rel.f_X_+=(2.*(exp(-pair.X_Mm)/4.+exp(-pair.X_MM)-P2) -var)/denom*count;
		rel.f_Y_+=(2.*(exp(-pair.Y_Mm)/4.+exp(-pair.Y_MM)-P2) -var)/denom*count;
	};
}

void 
inc_theta(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	if(P!=0){
		rel.theta_XY_+=( (exp(-pair.X_Mm-pair.Y_Mm)/4.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,2) )/pow(P*(1-P),2)*count;
	};
}

void 
inc_gamma(Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{
	float_t P=1-pair.m;
	if(P!=0 && P!=0.5){
		rel.gamma_XY_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/2.+exp(-pair.X_Mm-pair.Y_MM)/4.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1-P)*(rel.theta_XY_*(1+2*P)+rel.f_X_*P+P/(1-P) ) )/(P*(1-P) )/(0.5-P)/2*count;
		rel.gamma_YX_+=( 2.*(exp(-pair.X_Mm-pair.Y_Mm)/8.+exp(-pair.X_MM-pair.Y_Mm)/4.+exp(-pair.X_Mm-pair.Y_MM)/2.+exp(-pair.X_MM-pair.Y_MM) )-pow(P,3)-P*(1.-P)*(rel.theta_XY_*(1.+2.*P)+rel.f_Y_*P+P/(1.-P) ) )/(P*(1.-P) )/(0.5-P)/2.*count;
	};
}

void 
inc_Delta (Relatedness &rel, const Genotype_pair &pair, const size_t &count)
{

}

//TODO FIX THIS TERRIBLE HACK JOB YOU SCHMUCK!

//int eval;

bool failbound(const double &x)
{
    return (x<-1. || x>1.);
}

bool passbound(const double &x)
{
    return (x>=0  && x<=1.);
}

float_t
get_ll (const Relatedness &rel, const Genotype_pair &pair, const float_t count) 
{
//	eval++;
	/*This is the basic likelihood model that we are fitting. It essentially calculates the correlation coefficents
	for the first four moments of the joint distribution. The math needs to be cleaned up here. For now it is typed
	up to minimize the chance of typos. a, b, c and d are r.v. Representing the presence (1) or absence of (0) of
	the Major allele in the haploid genomes of individuals A (for a and b) and C (c and d). A is the r.v. defined as
	A~(a+b)/2 and C~(c+d)/2. Generally what we are doing here is calculating the expectations (e.g. E_A2) of the 
	joint distributions and using that to calculate the joint distribution itself. Problems can occur when 
	correlation coefficients are less than zero, because the probabilities of various observations (e.g. mm1mm2) can 
	become negative. These probabilities are forced to be zero, which may be a little arbitrary, but it seems to 
	work.*/

    if (failbound(rel.f_X_) || failbound(rel.f_Y_) || failbound(rel.theta_XY_) ) //|| failbound(rel.gamma_XY_) || failbound(rel.gamma_YX_) || failbound(rel.delta_XY_) || failbound(rel.Delta_XY_) )
        return nan("");

    if (!(passbound(pair.X_mm) && passbound(pair.Y_mm) && passbound(pair.X_Mm) && passbound(pair.Y_Mm) && passbound(pair.X_MM) && passbound(pair.Y_MM) ) ) {
        std::cerr << "Bad genotype...\n";
        return 0;//nan("");
    }

	real_t P, mm1mm2, Mm1mm2, MM1mm2, mm1Mm2, Mm1Mm2, MM1Mm2, mm1MM2, Mm1MM2, MM1MM2;

	P=1-pair.m;

	if (P==0 || P==1) {
        std::cerr << "Bad error...\n";
        return 0;
    }
	real_t e=0;//global_relatedness.e_Y_[freqtoi(pair.m)]-global_relatedness.e_X_[freqtoi(pair.m)];
	//P+=(global_relatedness.e_Y_[freqtoi(pair.m)]+global_relatedness.e_X_[freqtoi(pair.m)]);

	/*This comes from the inverse matrix of the one used to calculate the moments.*/
	
	real_t A=P+e;//+e(P)*P; 	//mean major allele frequency in A
	real_t C=P-e;//-e(P)*P; 	//mean major allele frequency in C
//	std::cerr << freqtoi(pair.m) << ", " << e << ", " << P << ", " << A << ", " << C << std::endl;	
	real_t Va=A*(1.-A);	//variance of the two haploid genomes of A 
	real_t Vc=C*(1.-C);	// "   "	"	" 	"     of B
	real_t Sa=sqrt(Va);	//standard deviation of haploid genomes A
	real_t Sc=sqrt(Vc);	// and "	"	"	"	C.

//`    std::cerr << Sa << ", " << Sc << std::endl;

	real_t E_A2  =(rel.f_X_*Va+2.*pow(A,2.)+Va)/2.; //Expectation of A^2
	real_t E_C2  =(rel.f_Y_*Vc+2.*pow(C,2.)+Vc)/2.; //
	real_t E_AC  =rel.theta_XY_*Sa*Sc+A*C;
//	if (Sa==0 || Sc==0) return 0;
	real_t ga=(1.-2.*A)/Sa;
	real_t gc=(1.-2.*C)/Sc;
	real_t E_A2C =(rel.gamma_XY_*Va*Sc*ga+A*A*C+Va*(rel.theta_XY_*(1.+2.*A)+rel.f_X_*C+C/(1-C) ) )/2.;
	real_t E_AC2 =(rel.gamma_YX_*Vc*Sa*gc+C*C*A+Vc*(rel.theta_XY_*(1.+2.*C)+rel.f_Y_*A+A/(1-A) ) )/2.;
	real_t ka=1./(1.-A)+1./A-3.;
	real_t kc=1./(1.-C)+1./C-3.;
	real_t E_A2C2=(rel.delta_XY_*sqrt(ka*kc)+rel.Delta_XY_)*Va*Vc+A*A*C*C+rel.f_X_*Va*C*C+rel.f_Y_*Vc*A*A+4.*rel.theta_XY_*Sa*Sc*A*A+C*2.*rel.gamma_XY_*Va*Sc*ga+2.*A*rel.gamma_YX_*Vc*Sa*gc;

	/*This comes from the inverse matrix of the one used to calculate the moments.*/

	mm1mm2=1.-6.*P+0.0*e+2*E_A2+2.*E_C2+8.0*E_AC-4.*E_A2C-4.*E_AC2+1.*E_A2C2;
	Mm1mm2=0.+4.*P+2.0*e-4*E_A2+0.*E_C2-10.*E_AC+8.*E_A2C+4.*E_AC2-2.*E_A2C2;
	MM1mm2=0.-1.*P-0.5*e+2*E_A2+0.*E_C2+2.0*E_AC-4.*E_A2C-0.*E_AC2+1.*E_A2C2;
	mm1Mm2=0.+4.*P-2.0*e+0*E_A2-4.*E_C2-10.*E_AC+4.*E_A2C+8.*E_AC2-2.*E_A2C2;
	Mm1Mm2=0.+0.*P+0.0*e+0*E_A2+0.*E_C2+12.*E_AC-8.*E_A2C-8.*E_AC2+4.*E_A2C2;
	MM1Mm2=0.+0.*P+0.0*e+0*E_A2+0.*E_C2-2.0*E_AC+4.*E_A2C+0.*E_AC2-2.*E_A2C2;
	mm1MM2=0.-1.*P+0.5*e+0*E_A2+2.*E_C2+2.0*E_AC+0.*E_A2C-4.*E_AC2+1.*E_A2C2;
	Mm1MM2=0.+0.*P+0.0*e+0*E_A2+0.*E_C2-2.0*E_AC+0.*E_A2C+4.*E_AC2-2.*E_A2C2;
	MM1MM2=0.+0.*P+0.0*e+0*E_A2+0.*E_C2+0.0*E_AC+0.*E_A2C+0.*E_AC2+1.*E_A2C2;

/*    if(P>0.98){

    std::cerr << "Par:" << rel.f_X_ << ", " << rel.f_Y_ << ", " << rel.theta_XY_ << ", " << rel.gamma_XY_ << ", " << rel.gamma_YX_ << ", " << rel.delta_XY_ << ", " << rel.Delta_XY_ << "\n";

    std::cerr << ga << ", " << ka << "\n";
    std::cerr << E_A2 << ", " << E_C2 << ", " << E_AC << ", " << E_A2C <<", " << E_AC2 << ", " << E_A2C2 << "\n";

    std::cerr << P << " : " << mm1mm2 << ", " << Mm1mm2 << ", " << MM1mm2 << "\n";
    std::cerr << P << " : " << mm1Mm2 << ", " << Mm1Mm2 << ", " << MM1Mm2 << "\n";
    std::cerr << P << " : " << mm1MM2 << ", " << Mm1MM2 << ", " << MM1MM2 << "\n";

    }*/

    /*
#define MIN -1.0e-16

    if (mm1mm2<MIN) return nan(""); 
	if (Mm1mm2<MIN) return nan(""); 
	if (MM1mm2<MIN) return nan(""); 
	if (mm1Mm2<MIN) return nan(""); 
	if (Mm1Mm2<MIN) return nan(""); 
	if (MM1Mm2<MIN) return nan(""); 
	if (mm1MM2<MIN) return nan(""); 
	if (Mm1MM2<MIN) return nan(""); 
	if (MM1MM2<MIN) return nan(""); 
    */
    /*	float_t S=pow(mm1mm2+Mm1mm2+MM1mm2+mm1Mm2+Mm1Mm2+MM1Mm2+mm1MM2+Mm1MM2+MM1MM2, 2);

	if(S>1){
		mm1mm2/=S;
		Mm1mm2/=S;
		MM1mm2/=S;
		mm1Mm2/=S;
		Mm1Mm2/=S;
		MM1Mm2/=S;
		mm1MM2/=S;
		Mm1MM2/=S;
		MM1MM2/=S;
	};
	
	float_t E[9];
*/
/*
        if (mm1mm2>0) E[0]=log(mm1mm2*pair.X_mm*pair.Y_mm);
     	else E[0]=-FLT_MAX;
	if (Mm1Mm2>0) E[1]=log(Mm1Mm2*pair.X_Mm*pair.Y_Mm);
     	else E[1]=-FLT_MAX;
        if (MM1MM2>0) E[2]=log(MM1MM2*pair.X_MM*pair.Y_MM);
     	else E[2]=-FLT_MAX;
        if (mm1Mm2>0) E[3]=log(mm1Mm2*pair.X_mm*pair.Y_Mm);
     	else E[3]=-FLT_MAX;
        if (MM1Mm2>0) E[4]=log(MM1Mm2*pair.X_MM*pair.Y_Mm);
     	else E[4]=-FLT_MAX;
        if (Mm1mm2>0) E[5]=log(Mm1mm2*pair.X_Mm*pair.Y_mm);
     	else E[5]=-FLT_MAX;
        if (Mm1MM2>0) E[6]=log(Mm1MM2*pair.X_Mm*pair.Y_MM);
     	else E[6]=-FLT_MAX;
        if (mm1MM2>0) E[7]=log(mm1MM2*pair.X_mm*pair.Y_MM);
     	else E[7]=-FLT_MAX;
        if (MM1mm2>0) E[8]=log(MM1mm2*pair.X_MM*pair.Y_mm);
     	else E[8]=-FLT_MAX;
	    std::sort(E, E+9);
*/
    float K=mm1mm2*pair.X_mm*pair.Y_mm+Mm1Mm2*pair.X_Mm*pair.Y_Mm+MM1MM2*pair.X_MM*pair.Y_MM+mm1Mm2*pair.X_mm*pair.Y_Mm+MM1Mm2*pair.X_MM*pair.Y_Mm+Mm1mm2*pair.X_Mm*pair.Y_mm+Mm1MM2*pair.X_Mm*pair.Y_MM+mm1MM2*pair.X_mm*pair.Y_MM+MM1mm2*pair.X_MM*pair.Y_mm;
    if (K<=0) {
//        std::cerr << "bad bounds: " << K << " " << count << std::endl;
        return nan("");
    }

/*    std::cerr << " - " << mm1mm2 << ", " << Mm1mm2 << ", " << MM1mm2 << "\n";
    std::cerr << " - " << mm1Mm2 << ", " << Mm1Mm2 << ", " << MM1Mm2 << "\n";
    std::cerr << " - " << mm1MM2 << ", " << Mm1MM2 << ", " << MM1MM2 << "\n";

    std::cerr << "bounds: " << mm1mm2 << ", " << P << ", " << e << ", " << E_A2 << ", " << E_C2 << ", " << E_AC << ", " << E_A2C<< K << " " << count << std::endl;
*/  return log(K)*count;


	//std::cerr << E[7] << ":" << E[8] << ":" << S  << ", " << count << std::endl;
	//return (log(1+exp(E[0]-E[8])+exp(E[1]-E[8])+exp(E[2]-E[8])+exp(E[3]-E[8])+exp(E[4]-E[8])+exp(E[5]-E[8])+exp(E[6]-E[8])+exp(E[7]-E[8]) )+E[8] )*count;
}


double
rel_ll (const gsl_vector *v, void *void_hashed_genotypes_p)
{
	Relatedness rel;
	std::vector < std::pair<Genotype_pair_tuple, size_t> > *hashed_genotypes_p=(std::vector <std::pair <Genotype_pair_tuple, size_t> > *) void_hashed_genotypes_p;
	Genotype_pair first;
	size_t count;

	rel.f_X_ = gsl_vector_get(v, 0);
	rel.f_Y_ = gsl_vector_get(v, 1);
	rel.theta_XY_ = gsl_vector_get(v, 2);
	rel.gamma_XY_ = gsl_vector_get(v, 3);
	rel.gamma_YX_ = gsl_vector_get(v, 4);
	rel.Delta_XY_ = gsl_vector_get(v, 5);
	rel.delta_XY_ = gsl_vector_get(v, 6);
	
	float_t sum=0;
	
	std::pair<Genotype_pair_tuple, size_t> *pair;
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator end=hashed_genotypes_p->end();
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator it=hashed_genotypes_p->begin();

	while(it!=end){
		first=Genotype_pair::from_tuple(it->first);
		count=it->second;
		if (first.m>0) sum+=get_ll(rel, first, count);
		it++;
	}
	if (std::isnan(sum) ) return nan("");
	return double(-sum);
}

void
rel_dll (const gsl_vector *v, void *void_hashed_genotypes_p,
               gsl_vector *df)
{
    gsl_vector_set_all(df, 0);
	std::vector < std::pair<Genotype_pair_tuple, size_t> > *hashed_genotypes_p=(std::vector <std::pair <Genotype_pair_tuple, size_t> > *) void_hashed_genotypes_p;
	Genotype_pair first;
	std::pair<Genotype_pair_tuple, size_t> *pair;
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator end=hashed_genotypes_p->end();
	std::vector<std::pair<Genotype_pair_tuple, size_t> >::iterator it=hashed_genotypes_p->begin();
	Constants <float_t, const std::pair<const Genotype_pair &, const Relatedness &> > consts(REL_CNTS, REL_ARRAY);
    double last_m=0;
	size_t count;
    Relatedness rel;

    rel.f_X_ = gsl_vector_get(v, 0);
	rel.f_Y_ = gsl_vector_get(v, 1);
	rel.theta_XY_ = gsl_vector_get(v, 2);
	rel.gamma_XY_ = gsl_vector_get(v, 3);
	rel.gamma_YX_ = gsl_vector_get(v, 4);
	rel.Delta_XY_ = gsl_vector_get(v, 5);
	rel.delta_XY_ = gsl_vector_get(v, 6);

    while (it!=end ){
		first=Genotype_pair::from_tuple(it->first);
	    if (first.m!=last_m){
		    consts.recalculate( std::pair <Genotype_pair, Relatedness> (first, rel) );
			last_m=first.m;
		}
		count=it->second;
        /*
        if (failbound(rel.f_X_) || failbound(rel.f_Y_) || failbound(rel.theta_XY_) || failbound(rel.gamma_XY_) || failbound(rel.gamma_YX_) || failbound(rel.delta_XY_) || failbound(rel.Delta_XY_) ){
			gsl_vector_set(df, 0, FLT_MAX );
			gsl_vector_set(df, 1, FLT_MAX );
			gsl_vector_set(df, 2, FLT_MAX );
			gsl_vector_set(df, 3, FLT_MAX );
			gsl_vector_set(df, 4, FLT_MAX );
			gsl_vector_set(df, 5, FLT_MAX );
			gsl_vector_set(df, 6, FLT_MAX );
            return;
        }*/
		if (first.m>0 && first.m<1.) {
			gsl_vector_set(df, 0, gsl_vector_get(df, 0)-H0(first, consts)*count );
			gsl_vector_set(df, 1, gsl_vector_get(df, 1)-H1(first, consts)*count );
			gsl_vector_set(df, 2, gsl_vector_get(df, 2)-H2(first, consts)*count );
			gsl_vector_set(df, 3, gsl_vector_get(df, 3)-H3(first, consts)*count );
			gsl_vector_set(df, 4, gsl_vector_get(df, 4)-H4(first, consts)*count );
			gsl_vector_set(df, 5, gsl_vector_get(df, 5)-H5(first, consts)*count );
			gsl_vector_set(df, 6, gsl_vector_get(df, 6)-H6(first, consts)*count );
        }
		++it;
	}
}

void
rel_fdf (const gsl_vector *x, void *params,
                double *f, gsl_vector *df)
{
      *f = rel_ll(x, params);
        rel_dll(x, params, df);
}

/*Maximizes the relatedness*/
bool 
maximize_gsl(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes, const int MAXITER)
{

	const gsl_multimin_fdfminimizer_type *T=gsl_multimin_fdfminimizer_conjugate_fr;
	//const gsl_multimin_fdfminimizer_type *T=gsl_multimin_fdfminimizer_conjugate_pr;
	//const gsl_multimin_fdfminimizer_type *T=gsl_multimin_fdfminimizer_vector_bfgs2; 
//	const gsl_multimin_fdfminimizer_type *T=gsl_multimin_fdfminimizer_vector_bfgs; 
//    const gsl_multimin_fdfminimizer_type *T=gsl_multimin_fdfminimizer_steepest_descent;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point and stepsizes. I would really prefer to do this with a 
	 * Newton-Raphson method, which just inexplicably like more than the 
	 * Nelder-Mead, but I'm being lazy today, or more accurately there are 
	 * fundamental problems with setting the problem up to use the NR method.
	 */
	gsl_vector *x=gsl_vector_alloc(7);
	gsl_vector *last_x=gsl_vector_alloc(7);

	gsl_vector_set(x, 0, rel.f_X_);
	gsl_vector_set(x, 1, rel.f_Y_);
	gsl_vector_set(x, 2, rel.theta_XY_);
	gsl_vector_set(x, 3, rel.gamma_XY_);
	gsl_vector_set(x, 4, rel.gamma_YX_);
	gsl_vector_set(x, 5, rel.Delta_XY_);
	gsl_vector_set(x, 6, rel.delta_XY_);

	gsl_vector *ss=gsl_vector_alloc(7);
	gsl_vector_set_all(ss,0.05);	
	gsl_multimin_function_fdf gsl_func;


	gsl_func.n=7;
	gsl_func.f = &rel_ll;
	gsl_func.df = &rel_dll;
	gsl_func.fdf = &rel_fdf; 
//
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );
    gsl_func.params=&hashed_genotypes_vector;

    /*
    struct multimin_params optim_par = {.1,1e-2,100,1e-3,1e-5,2,0};
    double minimum;
    double xmin[7] = {-1};
    double xmax[7] = {1};
    unsigned type[7]={3};
    multimin(7, x, &minimum, type, xmin, xmax, &rel_ll, &rel_dll, &rel_fdf, (void *) hashed_genotypes_vector, optim_par);
*/

    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 7);
    gsl_multimin_fdfminimizer *last_s = gsl_multimin_fdfminimizer_alloc (T, 7);
	gsl_multimin_fdfminimizer_set (s, &gsl_func, x, 0.001, 1e-4);//ss);
	
//	eval=0;
    
	do {
		iter++;
        gsl_vector_memcpy(last_x, s->x);
		status = gsl_multimin_fdfminimizer_iterate(s);
        
       for (int i = 0; i < 7; i++) printf ("%10.3f", gsl_vector_get (s->x, i));
        printf (" (s) f() = %7.3f size = %.3f\n", s->f, 7);
        for (int i = 0; i < 7; i++) printf ("%10.3f", gsl_vector_get (last_x, i));
        printf (" (l)...\n", last_s->f, 7);
        if(std::isnan(s->f) ){ //If this is nonsense go back and decrease step size...
            gsl_vector_memcpy(s->x, last_x);
   //         s->step_size=0.00001;
            gsl_multimin_fdfminimizer_restart(s);
        }

		if (status) break;

        for (int i = 0; i < 7; i++) printf ("%10.3f", gsl_vector_get (s->x, i));
        printf ("(a) f() = %7.3f size = %.3f\n", s->f, 7);
        //Test to see if step takes you to a good place, undo the step if it does not...

		//size = gsl_multimin_fdfminimizer_size (s);
		status = gsl_multimin_test_gradient (s->gradient, 1e-3);
	}  while (status == GSL_CONTINUE && iter < MAXITER);

    if (status == GSL_CONTINUE)
        std::cerr << "Maximum iterations exceded . . . " << iter << std::endl;
    rel.success_= (status != GSL_CONTINUE);
	rel.f_X_ = gsl_vector_get(s->x, 0);
	rel.f_Y_ = gsl_vector_get(s->x, 1);
	rel.theta_XY_ = gsl_vector_get(s->x, 2);
	rel.gamma_XY_ = gsl_vector_get(s->x, 3);
	rel.gamma_YX_ = gsl_vector_get(s->x, 4);
	rel.Delta_XY_ = gsl_vector_get(s->x, 5);
	rel.delta_XY_ = gsl_vector_get(s->x, 6);
	
	rel.max_ll_=rel_ll(x, &hashed_genotypes_vector);

    gsl_vector_free(x);
	gsl_vector_free(last_x);
    return (status != GSL_CONTINUE);
}

class small_rel{

	public:
	size_t skip;
	small_rel(size_t skip_) {
		skip=skip_;
	}
	double
	rel_ll_2 (const gsl_vector *v, void *void_hashed_genotypes_p)
	{
       		gsl_vector *x=gsl_vector_alloc(7);

		for (size_t y=0; y<7;y++){
			if(y==skip){
   		     		gsl_vector_set(x, y, 0);
			} else {
				gsl_vector_set(x, y, gsl_vector_get(v, y-(y<=skip) ) );
			}
		}
		return rel_ll(v, void_hashed_genotypes_p);
	}
};

void
get_llr(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> hashed_genotypes)
{
	const gsl_multimin_fminimizer_type *T=gsl_multimin_fminimizer_nmsimplex;
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *ss, *x;
//..	gsl_vector_set_all(ss,0.15);	
	gsl_multimin_function gsl_func;
	
	size_t iter = 0;
	int status;
	double size;

	std::vector <float_t> values={rel.f_X_, rel.f_Y_, rel.theta_XY_,rel.gamma_XY_,rel.gamma_YX_, rel.Delta_XY_, rel.delta_XY_};
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );


	/* Starting point and stepsizes. I would really prefer to do this with a 
	 * Newton-Raphson method, which just inexplicably like more than the 
	 * Nelder-Mead, but I'm being lazy today, or more accurately there are 
	 * fundamental problems with setting the problem up to use the NR method.
	 */
	for (size_t w=0; w<7; ++w) {

//		small_rel this_rel(w);
		
//		ss=gsl_vector_alloc(7);
		x=gsl_vector_alloc(7);
	    gsl_vector_set_all(x,0.0);	

		rel.null_ll_ = rel_ll(x, &hashed_genotypes_vector);

	        gsl_vector_set(x, 0, rel.f_X_);
	        gsl_vector_set(x, 1, rel.f_Y_);
	        gsl_vector_set(x, 2, rel.theta_XY_);
	        gsl_vector_set(x, 3, rel.gamma_XY_);
	        gsl_vector_set(x, 4, rel.gamma_YX_);
	        gsl_vector_set(x, 5, rel.Delta_XY_);
	        gsl_vector_set(x, 6, rel.delta_XY_);

		rel.max_ll_ = rel_ll(x, &hashed_genotypes_vector);



//TODO FIX ME!
/*		gsl_func.n=6;
		gsl_func.f = rel_ll;
		gsl_func.params=&hashed_genotypes;

		s = gsl_multimin_fminimizer_alloc (T, 6);
		gsl_multimin_fminimizer_set (s, &gsl_func, x, ss);
		do {
			iter++;
			status = gsl_multimin_fminimizer_iterate(s);
      
			if (status) break;

			size = gsl_multimin_fminimizer_size (s);
			status = gsl_multimin_test_size (size, 1e-4);

		}  while (status == GSL_CONTINUE && iter < 600);
*/

		switch (w) {
			case 0:
	        		gsl_vector_set(x, 0, 0);
				rel.f_X_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.f_X_);
			break;
			case 1:
	        		gsl_vector_set(x, 1, 0);
				rel.f_Y_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.f_Y_);
			break;
			case 2:
	        		gsl_vector_set(x, 2, 0);
				rel.theta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.theta_XY_);
			break;
			case 3:
	        		gsl_vector_set(x, 3, 0);
				rel.gamma_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.gamma_XY_);
			break;
			case 4:
	        		gsl_vector_set(x, 4, 0);
				rel.gamma_YX_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.gamma_YX_);
			break;
			case 5:
	        		gsl_vector_set(x, 5, 0);
				rel.Delta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.Delta_XY_);
			break;
			case 6:
	        		gsl_vector_set(x, 6, 0);
				rel.delta_XY_ll = rel_ll(x, &hashed_genotypes_vector)-rel.max_ll_;
	        		gsl_vector_set(x, 0, rel.delta_XY_);
			break;
		}
	}
	
}

void
get_95CI(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> hashed_genotypes)
{	
    std::vector <float_t> values={rel.f_X_, rel.f_Y_, rel.theta_XY_,rel.gamma_XY_,rel.gamma_YX_, rel.Delta_XY_, rel.delta_XY_};
	std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );

    gsl_matrix *invJ= gsl_matrix_alloc(7,7);
    get_invJ(invJ, rel, hashed_genotypes);

    rel.f_X_ll=1.96*sqrt(gsl_matrix_get(invJ,0,0) );	
    rel.f_Y_ll=1.96*sqrt(gsl_matrix_get(invJ,1,1) );	
    rel.theta_XY_ll=1.96*sqrt(gsl_matrix_get(invJ,2,2) );	
    rel.gamma_XY_ll=1.96*sqrt(gsl_matrix_get(invJ,3,3) );	
    rel.gamma_YX_ll=1.96*sqrt(gsl_matrix_get(invJ,4,4) );	
    rel.Delta_XY_ll=1.96*sqrt(gsl_matrix_get(invJ,5,5) );	
    rel.delta_XY_ll=1.96*sqrt(gsl_matrix_get(invJ,6,6) );	
		
    gsl_vector *x=gsl_vector_alloc(7);
    gsl_vector_set_all(x,0.0);	

	rel.null_ll_ = rel_ll(x, &hashed_genotypes_vector);

    gsl_vector_set(x, 0, rel.f_X_);
    gsl_vector_set(x, 1, rel.f_Y_);
    gsl_vector_set(x, 2, rel.theta_XY_);
    gsl_vector_set(x, 3, rel.gamma_XY_);
    gsl_vector_set(x, 4, rel.gamma_YX_);
    gsl_vector_set(x, 5, rel.Delta_XY_);
    gsl_vector_set(x, 6, rel.delta_XY_);

	rel.max_ll_ = rel_ll(x, &hashed_genotypes_vector);
}

#ifdef NOMPI
int estimateRel(int argc, char *argv[])
{
	std::cerr << "Warning: this program should generate AICc's. However, it doesn't and its _ll variables don't mean that much. Also it is way too slow.\n"; 

	/* All the variables that can be set from the command line */

	std::string gcf_name="", rel_name="";

	bool l2o=false, called=false, use_newton=false;
    int startx=0, starty=1;

	Environment env;
	env.set_name("mapgd relatedness");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Uses a maximum likelihood approach to estimate pairwise relatedness.");

	env.optional_arg('i', "input", 	gcf_name, "an error occurred while displaying the help message.", "input filename (default stdout)");
	env.optional_arg('o', "output", rel_name, "an error occurred while displaying the help message.", "output filename (default stdin)");
	env.optional_arg('x', "x", startx, "an error occurred while displaying the help message.", "fist individual x (default 0)");
	env.optional_arg('y', "y", starty, "an error occurred while displaying the help message.", "fist individual y (default 1)");
	env.flag(	'n', "newton", 	&use_newton, 		&flag_set, "an error occurred while displaying the help message.", "uses the NR optimization. Doesn't work.");
	env.flag(	'l', "l2o", 	&l2o, 		&flag_set, "an error occurred while displaying the help message.", "uses the 'leave 2 out' procedure of calculating allele freq.");
	env.flag(	'c', "called", 	&called,	&flag_set, "an error occurred while displaying the help message.", "ignore genotypic likelihoods and treat data as error free.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the help message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Population> gcf_in; 	// Open the file with genotypic probabilities.
	Indexed_file <Population> gcf_mem; 	// the in memory gcf_file.
	std::stringstream file_buffer;

	Flat_file <Relatedness> rel_out; 		// Open a file for output.

	Relatedness relatedness;		//The class to read to
	Population genotype;			//The class to write from

	if (gcf_name.size()!=0)
		gcf_in.open(gcf_name.c_str(), READ );
	else
		gcf_in.open( READ );

	if (rel_name.size()!=0)
		rel_out.open(rel_name.c_str(), WRITE );
	else
		rel_out.open( WRITE );
	
	genotype=gcf_in.read_header();			//This gives us the sample names.

	gcf_mem.open(&file_buffer, BINARY | WRITE );
	gcf_mem.set_index(gcf_in.get_index() );
	gcf_mem.write_header(genotype);

	while(gcf_in.table_is_open() ){
		gcf_in.read(genotype);
		gcf_mem.write(genotype);
	}
	gcf_mem.close();

	rel_out.write_header(relatedness);

    std::vector <std::string> name=genotype.get_sample_names();

	std::map <Genotype_pair_tuple, size_t> hashed_genotypes;
	std::map <Genotype_pair_tuple, size_t> down_genotypes;
	
	size_t sample_size=genotype.get_sample_names().size();
	

	for (size_t x=startx; x<sample_size; ++x){
		for (size_t y=x+1; y<sample_size; ++y){
            if (x==startx and y < starty) y=starty;
			relatedness.set_X_name(x);
			relatedness.set_Y_name(y);
			hashed_genotypes=hash_genotypes(file_buffer, x, y, l2o, called);
/*
            gsl_vector *t=gsl_vector_alloc(7);

	gsl_vector_set(t, 0, 
            0.607);
	gsl_vector_set(t, 1, 
            0.607);
	gsl_vector_set(t, 2, 
            0.0); 
	gsl_vector_set(t, 3, 
            0.0);
	gsl_vector_set(t, 4, 
            0.0);
	gsl_vector_set(t, 6, 
            0.0);
	gsl_vector_set(t, 5, 
            0.369);
            std::vector <std::pair <Genotype_pair_tuple, size_t> > hashed_genotypes_vector(hashed_genotypes.begin(), hashed_genotypes.end() );
            std::cerr << rel_ll(t, &hashed_genotypes_vector);
            return 0;
*/

            //			down_genotypes=downsample_genotypes(file_buffer, x, y);
//			down_genotypes=hash_genotypes(file_buffer, x, y, l2o);
			relatedness.zero();
			//set_e(relatedness, hashed_genotypes);
			gestimate(relatedness, hashed_genotypes);
#ifdef EIGEN
			//newton(relatedness, down_genotypes);
			maximize_gsl(relatedness, down_genotypes, ITER_MAX);
#else
//			if( !newton(relatedness, hashed_genotypes) ) relatedness.success_= true;
		if( maximize_gsl(relatedness, hashed_genotypes, ITER_MAX) ) relatedness.success_= true;
            	else relatedness.success_= false;
#endif
            if (relatedness.success_) get_95CI(relatedness, hashed_genotypes);
			relatedness.set_X_name(name[x]);
			relatedness.set_Y_name(name[y]);
			rel_out.write(relatedness);
        }
	}
	rel_out.close();
	return 0;					//Since everything worked, return 0!.
}
#endif

#else 

int 
estimateRel(int argc, char *argv[])
{
	std::cerr << "This command depends on gsl, which could not be found. Please e-mail matthew.s.ackerman@gmail.com for help.\n";
	return 0;
}

#endif 
