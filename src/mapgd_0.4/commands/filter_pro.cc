/* 

command filter:

*/

#include "filter_pro.h"

double
get_coverage_z_score(const Locus &locus, const std::vector <double> &rates)
{
	double sum=0, norm=0;

	std::vector <quartet_t>::const_iterator q_it=locus.sample.cbegin();
	std::vector <quartet_t>::const_iterator q_end=locus.sample.cend();

	std::vector <double>::const_iterator rate_it=rates.cbegin();
    std::cerr << locus.get_abs_pos();
	while (q_it!=q_end)
	{
        std::cerr << ", " << count(*q_it); 
		if (*rate_it > 0)
		{

            sum += pow(gsl_cdf_ugaussian_Pinv(gsl_cdf_poisson_P(count(*q_it), *rate_it) ), 2.0);		//This should be the probability of sampling count() reads from a Poisson distribution with rate *rate_it
			norm += 1;//0.25;	//This is the mean of the pdf with rate *rate_it
		}
		q_it++;
		rate_it++;
	}
	std::cerr << ", " << sum << ", " << norm << std::endl; 
	return 1.-gsl_cdf_chisq_P(sum, norm);
}

double 
modifided_negative_binomial_pdf(const int &k, const double &p, const double &r )
{
    return powf(p,k)*powf(1.-p,r)*gsl_sf_gamma(k+r)/gsl_sf_gamma(k+1)/gsl_sf_gamma(r);
}

double
modifided_cdf_negative_binomial_P(const int &k, const double &p, const double &r )
{
    return 1.-gsl_sf_beta_inc(r,k,p);
}

double
get_coverage_z_score_neg_bin(const Locus &locus, const std::vector < std::array < double, 2 > > &par)
{
	double sum=0, norm=0;

	std::vector <quartet_t>::const_iterator q_it=locus.sample.cbegin();
	std::vector <quartet_t>::const_iterator q_end=locus.sample.cend();

	std::vector < std::array<double, 2> >::const_iterator p_it=par.cbegin();

	while (q_it!=q_end)
	{
		if ( (*p_it)[0] > 0 & (*p_it)[0] < 1 )
		{
            double p=(*p_it)[0], r=(*p_it)[1];
			sum += pow(gsl_cdf_ugaussian_Pinv(modifided_cdf_negative_binomial_P(count(*q_it), (*p_it)[0] , (*p_it)[1] ) ), 2.0);		//This should be the probability of sampling count() reads from a Poisson distribution with rate *rate_it
            std::cerr << locus.get_abs_pos() << ", " << (*p_it)[0] << ", " << (*p_it)[1] << ", E: " << p*r/(1-p) << ", O: " <<  count(*q_it) << ", " << modifided_negative_binomial_pdf(count(*q_it), (*p_it)[0] , (*p_it)[1] ) << ", " << modifided_cdf_negative_binomial_P(count(*q_it), (*p_it)[0] , (*p_it)[1] ) << std::endl;		//This should be the probability of sampling count() reads from a Poisson distribution with rate *rate_it
			norm += 1;//0.25;	//This is the mean of the pdf with rate *rate_it
		}
		q_it++;
		p_it++;
	}
//	std::cerr << sum << ", " << norm << std::endl; 
	return 1.-gsl_cdf_chisq_P(sum, norm);
}

std::vector< std::array <double, 2 > >
gc_adjust_lin_2(const std::vector <double> &old_rates, const double &gc, const double &x0, const double &x1, const double &b0, const double &b1)
{
	std::vector <std::array<double, 2> > rates;//(old_rates.size(),{0,0});
    for (int x=0; x<old_rates.size(); x++)
    {
        std::array<double, 2> pair { {old_rates[x], 0} };
        rates.push_back(pair);
    }

    //=old_rates;
	std::vector < std::array<double, 2 > >::iterator it=rates.begin();
	std::vector < std::array<double, 2 > >::iterator end=rates.end();

    double norm=320.344843;

	while (it!=end)
	{
		(*it)[1]=(*it)[0]*(x1+gc*b1)/norm;
		(*it)[0]=(x0+gc*b0);
		it++;
	}
	return rates;	
}

std::vector<double>
gc_adjust_lin(const std::vector <double> &old_rates, const double &gc, const double &a, const double &b)
{
	std::vector <double> rates=old_rates;
	std::vector <double>::iterator it=rates.begin();
	std::vector <double>::iterator end=rates.end();

	while (it!=end)
	{
		(*it)=(*it)*(a+gc*b);
		it++;
	}
	return rates;	
}

std::vector<double>
gc_adjust_norm(const std::vector <double> &old_rates, const double &gc, const double &mean, const double stdev)
{
	std::vector <double> rates=old_rates;
	std::vector <double>::iterator it=rates.begin();
	std::vector <double>::iterator end=rates.end();
	double var=powf(stdev,2);
	while (it!=end)
	{
		(*it)=1./sqrt(2*M_PI*var)*exp(-powf(gc-mean,2)/(2*var) )*(*it);
	}
	return rates;
}


int filter_pro(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::vector <double> rates;

	std::string in_file="", out_file="";
	
	float_t a=0, b=0, c=0, d=0, min_z=0;

	bool binary=false, linear=false, normal=false, fit=false;

	int max_coverage=CNT_MAX, min_coverage=4, w=300;

	Environment env;
	env.set_name("mapgd filterpro");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Filter sites in '.pro' files based on criteria.");
	
	env.optional_arg('i',"input", 	in_file, "please provide a filename.", "the name of the input file (default stdin).");
	env.optional_arg('o',"output", 	out_file, "please provide a filename.", "the name of an output file (default stdout).");

	env.optional_arg('c',"min-cov", min_coverage, "please provide an integer.", "minimum coverage for a population at a site for the site to be used (default 4).");
	env.optional_arg('C',"max-cov", max_coverage, "please provide an int.", "max coverage for a population at a site for the site to be used (default CNT_MAX).");


	env.optional_arg('a',"", a, "please provide an float.", "parameter a (intercept for linear, mean for normal, default none).");
	env.optional_arg('b',"", b, "please provide an float.", "parameter b (slope for linear, stdev for normal, default none).");
	env.optional_arg('d',"", c, "please provide an float.", "parameter a (intercept for linear, mean for normal, default none).");
	env.optional_arg('e',"", d, "please provide an float.", "parameter b (slope for linear, stdev for normal, default none).");


    env.optional_arg('w',"", w, "please provide an int.", "window size (used for linear and normal GC model. Default 300).");

	env.optional_arg('m',"rate-cov", rates, "please provide a float.", "list of rate coverages for the population (default none).");
	env.optional_arg('z',"min-z", min_z, "please provide an float.", "minimum p-value for accepting a site based on coverage model (default CNT_MAX).");

	env.flag(	'f',"fit", 	&fit,		&flag_set, 	"please provide an int.", "don't filter sites, instead estimate model parameters.");
	env.flag(	'l',"linear", 	&linear,	&flag_set, 	"please provide an int.", "apply linear GC coverage model.");
//	env.flag(	'2',"linear2", 	&linear,	&flag_set, 	"please provide an int.", "apply linear GC coverage model.");
	env.flag(	'n',"normal", 	&normal,	&flag_set, 	"please provide an int.", "apply normal GC coverage model.");

	env.flag(	'b',"binary", 	&binary,	&flag_set, 	"please provide an int.", "output in binary mode (fast).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Locus s, *buffer_locus, *s_half, *s_it, *s_end;
	double *gc, *g_it, *g_half, *g_end;

	Indexed_file <Locus> pro_in;
	Indexed_file <Locus> pro_out;

	if (in_file.size()==0)	pro_in.open(std::fstream::in);
	else pro_in.open(in_file.c_str(), std::fstream::in);

	if (out_file.size()==0)	pro_out.open(std::fstream::out);
	else pro_out.open(out_file.c_str(), std::fstream::out);

	s=pro_in.read_header();

	pro_out.set_index(pro_in.get_index() );
	pro_out.write_header(s);

	if (rates.size()==1)
	{
		double value=rates[0];
		rates.resize(s.sample.size() );
		std::fill(rates.begin(), rates.end(), value);
	}

	int GC=0, N=0;

	buffer_locus=new Locus[w];
	s_end=buffer_locus+w;
	std::fill_n(buffer_locus, w, s);

	s_it=buffer_locus;
	s_end=buffer_locus+w;

	while (s_it != s_end && pro_in.table_is_open() )
	{
		pro_in.read(*s_it);
		if (s_it->ref==Base::ctob('G') || s_it->ref==Base::ctob('C') ) 
			GC++;
		if (s_it->ref==Base::ctob('N')  ) 
			N++;
		s_it++;
	}

	gc=new double [w];
	if (N<w)
		std::fill_n(gc, w, (GC)/(w-N) );
	else
		std::fill_n(gc, w, 0 );

	g_end=gc+w;

	s_it=buffer_locus;
	s_half=buffer_locus+w/2;
	g_it=gc;
	g_half=gc+w/2;
	double SX=0, SXsq=0, SY=0, SYsq=0, SXY=0, n=0;
	std::vector <int> sums(s.sample.size(), 0);
	
	while( pro_in.table_is_open() )
	{
		//std::cerr << s_it->get_abs_pos() << ", " << s_it->getcoverage() << std::endl;
		//pro_out.write(*s_it);
		if ( s_half->getcoverage() >= min_coverage && s_half->getcoverage() <= max_coverage)	
		{
			if (fit)
			{
				double cov = s_half->getcoverage();
				double gc  = *g_half;
				SY   += cov;
				SYsq += cov*cov;
				SX   += gc;
				SXsq += gc*gc;
				SXY  += gc*cov;
				n += 1.;

				std::vector <quartet_t>::const_iterator q_it=s_half->sample.cbegin();
				std::vector <quartet_t>::const_iterator q_end=s_half->sample.cend();
				std::vector <int>::iterator sum_it=sums.begin();
				std::vector <int>::iterator sum_end=sums.end();

				while (q_it!=q_end)
				{
					*sum_it+=count(*q_it);
					q_it++;
					sum_it++;
				}

				if (s_half->get_abs_pos() % 1000 == 0){
					std::cerr << "a:" << (SY*SXsq-SX*SXY)/(n*SXsq-SX*SX)/(SY/n) << ", " << "b:" << (n*SXY-SX*SY)/(n*SXsq-SX*SX)/(SY/n) << std::endl;
					std::cerr << "p (a):" << (SY*SXsq-SX*SXY)/(n*SXsq-SX*SX)/(SY/n) << ", " << "p (b):" << (n*SXY-SX*SY)/(n*SXsq-SX*SX)/(SY/n) << std::endl;
					std::cerr << "r (a):" << (SY*SXsq-SX*SXY)/(n*SXsq-SX*SX)/(SY/n) << ", " << "r (b):" << (n*SXY-SX*SY)/(n*SXsq-SX*SX)/(SY/n) << std::endl;
					sum_it=sums.begin();
					while (sum_it!=sum_end)
					{
						std::cerr << double(*sum_it)/double(n) << ",";
						sum_it++;
					}
					std::cerr << std::endl;
				}
			} else {
				if (normal)
				{
					std::vector <double> ad_rates=rates;
					if (get_coverage_z_score(*s_half, gc_adjust_lin(ad_rates, *g_it, a, b) ) > min_z)
					{
						pro_out.write(*s_half);
					}
				}
				else if (linear)
				{
					//std::vector <double> ad_rates=gc_adjust_lin(rates, *g_it, a, b);
					//std::cerr << s_half->get_abs_pos() << ", " << get_coverage_z_score(*s_half, ad_rates)  << ", " << s_half->getcoverage() << ", " << *g_half << std::endl;
					std::vector < std::array <double, 2 > > ad_rates=gc_adjust_lin_2(rates, *g_it, a, b, c, d);
					std::cerr << s_half->get_abs_pos() << ", " << get_coverage_z_score_neg_bin(*s_half, ad_rates)  << ", " << s_half->getcoverage() << ", " << *g_half << std::endl;
					if (get_coverage_z_score_neg_bin(*s_half, ad_rates) > min_z)
					{
						pro_out.write(*s_half);
					}
				}
				else{
	//				std::cerr << s_half->get_abs_pos() << ", " << get_coverage_z_score(*s_half, rates) << ", " << s_half->getcoverage() << ", " << *g_half << std::endl;
					if (get_coverage_z_score(*s_half, rates) > min_z)
					{
						pro_out.write(*s_half);
					} 
				}
			}
		}

		if (s_it->ref==Base('G') || s_it->ref==Base('C') ) 
			GC--;
		if (s_it->ref==Base('N')  ) 
			N--;
		pro_in.read(*s_it);
		*g_half=double(GC)/double(w-N);
		if (s_it->ref==Base('G') || s_it->ref==Base('C') ) 
			GC++;
		if (s_it->ref==Base('N')  ) 
			N++;
		if (++s_it==s_end)
			s_it=buffer_locus;
		if (++g_it==g_end)
			g_it=gc;
		if (++g_half==g_end)
			g_half=gc;
		if (++s_half==s_end)
			s_half=buffer_locus;
	}
	pro_out.close();
	return 0;					//Since everything worked, return 0!.
}
