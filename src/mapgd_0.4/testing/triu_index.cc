#include "triu_index.h"

Triangular_index::Triangular_index(size_t size)
{
	n_=size;
	size_= (n_*(n_-1)/2);
	k_=0;
	x_=0;
	y_=1;
	dirty_=false;
}

void 
Triangular_index::get_xy(size_t &x, size_t &y) 
{
	if (dirty_){
		if (k_<=n_){
			x_ = n_ - 2 - floor(sqrt(-8*k_ + 4*n_*(n_-1)-7)/2.0 - 0.5);
			y_ = k_ + x_ + 1 - n_*(n_-1)/2 + (n_-x_)*((n_-x_)-1)/2;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Array index out of bounds.\n"), __FILE__, __LINE__);
		}
		dirty_=false;
	}
	x=x_;
	y=y_;
}

size_t 
Triangular_index::get_k(void) const
{
	return k_;
}

size_t 
Triangular_index::get_k(const size_t &x, const size_t &y) 
{
	set(x,y);
	return k_;
}

size_t 
Triangular_index::size(void) const
{
	return size_;
}

void 
Triangular_index::set(const size_t &x, const size_t &y)
{
	x_=x;
	y_=y;
	k_ = (n_*(n_-1)/2) - (n_-x)*((n_-x)-1)/2 + y - x - 1;
}

Triangular_index& Triangular_index::operator++()
{
	k_++;
	if (++y_ >= n_)
	{
		y_=(++x_)+1;
	}
	return *this;
}

Triangular_index Triangular_index::operator++(int)
{
	Triangular_index tmp(*this); // copy
	operator++(); // pre-increment
	return tmp;   // return old value
}
