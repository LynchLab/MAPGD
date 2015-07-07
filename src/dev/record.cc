//classes for population genomic data. This file contains all of the data we could possible want  

//@breif: The basic file like class. It handels reading and writing "ROW" data.
template <class T> class record {
private:
	index index_;						//handels encoding/decoding id0.
	uint16_t id0;						//.
	uint64_t id1;						//uint64_t so we can index <i> Paris japonica </i>
	T data_;

	inline void writeb (std::ostream&) const;
	inline void writet (std::ostream&) const;
	inline void readb (std::istream&);
	inline void readt (std::istream&);
public:
	T get(void) const;
	T set(const T &);

	std::string getid0(void) const;				//
	uint32_t getid1(void) const;				//
	//Accessors and Mutators
	record & operator = (const record  &);

	friend std::ostream& operator << (std::ostream&, const record &) const;
	friend std::istream& operator >> (std::ostream&, record &);

	bool operator < (record &lhs, record &rhs) const;
}

class allele {
private:
	uint8_t minor_;						//idenity of minor allele (member).
	uint8_t major_;						//idenity major allele (member).
	uint8_t ploidy_;					//ploidy of the ?
	float_t *allele_frequency_; 				//all

	inline void writeb (std::ostream&) const;
	inline void writet (std::ostream&) const;
	inline void readb (std::istream&);
	inline void readt (std::istream&);
public:
	//constructor
	allele(const uint8_t &);
	//Accessors and Mutators
	const uint8_t & minor(void) const;			//access idenity of minor allele.
	const uint8_t & set_minor(const float_t &);		//set idenity of minor allele.
	const uint8_t & major(void) const;			//access idenity major allele.
	const uint8_t & set_major(void) const;			//set idenity major allele.
	float_t & frequency(void) const;			//acess frequency of major allele.
	float_t & set_frequency(const float_t &);		//set frequency of major allele.

	record & operator=(const allele  &);

	inline void write (std::ostream&) const;
	inline void read (std::istream&) ;
}

class haploid : public allele {
	inline void writeb (std::ostream&) const;
	inline void writet (std::ostream&) const;
	inline void readb (std::istream&);
	inline void readt (std::istream&);
public:
	//constructor
	//Accessors and Mutators
	float_t & genotype_frequency(const uint32_t &);		//return the frequency of the minor (0) or Major (1) genotype.
	float_t & M(void);					//return the frequency of the (M)ajor genotype.
	float_t & m(void);					//return the frequency of the (m)inor genotype.
	float_t & set_genotype_frequency(const uint32_t &);	//set the frequency of the minor (0) or Major (1) genotype.
	//?
	float_t & heterozygosity(void) const;			//heterozygosity.
	//?
	haploid & operator=(const hapoid  &);
	//IO
	inline void write (std::ostream&) const;
	inline void read (std::istream&) ;
}

class diploid : public allele {
private:
public:
	float_t genotype_frequency(const uint32_t &) const;	//
	float_t genotype_frequency(const uint32_t &, const uint32_t &) const;

	float_t set_genotype_frequency(const uint32_t *);	//set the frequency of the minor (0) or Major (1) genotype.

	float_t MM(void) const;
	float_t Mm(void) const;
	float_t mm(void) const;

	float_t fixation_index(void) const;		//HW statistic.
	float_t heterozygosity(void) const;		//heterozygosity.
}	

template <class T> class likelihood_statistics 
{
private:
	T allele_;		//The allele the likelihood_statistics describe.
	float_t fit_;		//log likelihood value of the fit.
	float_t error_;		//ml_error rate.
public:
	//constructor.
	//Accessors.
	float_t fit();
	float_t error();
	inline void write (std::ostream&) const;
	inline void read (std::istream&) ;
}

{
	uint32_t excluded;		//
	uint32_t coverage;		//population coverage.
	uint32_t ?;			//number of individual at the site.
	float_t effective_chromosomes;	//number of 'effective' chromosomes in the sample.

template <class T> class individual 
{
private:
	//COLUMN
	std::string name_;
	//ROW
	profile *profile_;
	likelihood_statistics <T> genotypes_; 
public:
	inline void write (std::ostream&) const;
	inline void read (std::istream&) ;
}

template <class T> class population  
{ 
private:
	uint32_t size_;
	//COLUMN
	std::string name_;
	//ROW
	T allele_;
	likelihood_statistics <T> monomorphic_, <T> polymorphic_; //hardy_weinberg_
	individual *individuals_;
public:
	std::string get_name(void) const;
	void set_name(std::string);

	float_t get_allele_frequency(void);

	inline void write (std::ostream&) const;
	inline void read (std::istream&) ;
}

/* A meta_poppulation is a collection of populations. This is the largest structure which we can 
 * easily assign indexs in our id0+id1 naming scheme. */
class meta_population  			
{
private:
	uint32_t size_;
	//COLUMN
	std::string name_;
	//ROW
	population *populations_;
	likelihood_statistics <?> Fst_;
public:
	std::string get_name(void) const;
	void set_name(std::string);
	inline void write (std::ostream&) const;
	inline void read (std::istream&);

	float_t get_Fst(void);
	float_t get_Fst_delta(void);
}
