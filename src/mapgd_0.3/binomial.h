#ifndef BINOMIAL_H_
#define BINOMIAL_H_

class binomial {
private:
	std::vector <uint32_t> fact_vector;
	float_t p_;
	uint32_t *bp_;
	uint32_t size_;
public:	
	binomial () {};
	binomial (const float_t &p) {p_=p;};
	~binomial (void){ fact_vector.clear(); }	//creates a function that returns log probabilites from a multinomial distribution with parameters float_t . . .

	uint32_t fact(const uint32_t&) ;		//returns the probabiltiy of the multinomial distribution . . .
	uint32_t binomial_coef(const uint32_t&, const uint32_t&);	//returns the log factorial of the count type numbers in the array
};

#endif
