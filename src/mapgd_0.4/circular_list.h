#ifndef _CICULAR_LIST_H_
#define _CICULAR_LIST_H_

#include <cstddef>
#include <iterator>
#include <iostream>

/*This is an extreamly poor and unstable implementation, inteneded only as a learning exercise. Please do not use this class*/

template <class T>
class Circular_list;

template <class T>
class C_it : public std::iterator<std::bidirectional_iterator_tag, T> 
{
private:
	struct Node {
		T value;
		Node *next;
		Node *prev;
	};
	Node *itr_;
public:
	C_it ()
	{
		itr_=new Node;
	}
	C_it (const T &val)
	{
		itr_=new Node;
		itr_->value=val;
	}

	void
	swap(C_it& other)
	{
		using std::swap;
		swap(itr_, other.itr_);
	}

	C_it & 
	operator++ (void) // Pre-increment
	{
		std::cout << "pre inc\n";
		itr_ = itr_->next;
		std::cout << itr_->value << std::endl;
		return *this;
	}

	C_it 
	operator++ (int) // Post-increment
	{
		std::cout << "post inc\n";
		C_it tmp(*this);
		itr_ = itr_->next;
		return tmp; 
	}

	C_it & 
	operator-- (void) // Pre-increment
	{
		itr_ = itr_->prev;
		return *this;
	}

	C_it 
	operator-- (int) // Post-increment
	{
		C_it tmp(*this);
		itr_ = itr_->prev;
		return tmp;
	}

	// two-way comparison: v.begin() == v.cbegin() and vice versa
	template <class Other_type>
	bool 
	operator == (const C_it<Other_type>& rhs) const
	{
		return itr_ == rhs.itr_;
	}

	template <class Other_type>
	bool 
	operator != (const C_it<Other_type>& rhs) const
	{
		return itr_ != rhs.itr_;
	}

	T & 
	operator* () const
	{
		return (itr_->value);
	}

	Node & 
	operator-> () const
	{
		return itr_;
	}

	operator C_it <const T>() const
	{
		return C_it <const T> (itr_);
	}

	friend Circular_list<T>;
};

//Cannot insert into empty list?
template <class T>
class Circular_list{
private:
	C_it<T> origin_;
public:
	typedef C_it<T> iterator;
	typedef const C_it<T> const_iterator;
	Circular_list()
	{
	       	origin_.itr_->next=origin_.itr_;
	       	origin_.itr_->prev=origin_.itr_;
	}

	Circular_list(const T &val)
	{
		origin_.itr_->value=val;			
	       	origin_.itr_->next=origin_.itr_;
	       	origin_.itr_->prev=origin_.itr_;
	}

	iterator insert (iterator &it, const T &val)
	{
		std::cout << "inserting " << val << std::endl;
        	C_it<T> *new_it = new C_it<T>(val);
		std::cout << **new_it << std::endl;
		new_it->itr_->prev=it.itr_->prev;
	       	new_it->itr_->next=it.itr_;
		it.itr_->prev->next=new_it->itr_;
		it.itr_->prev=new_it->itr_;
		std::cout << "prev is now " << it.itr_->prev->value << std::endl;
		return *new_it;
	}

	iterator insert (iterator &it, const size_t &n, const T &val)
	{
		std::cout << "inserting " << val << std::endl;
        	C_it<T> *new_it = new C_it<T>[n];
        	C_it<T> *first_it = new_it;
        	C_it<T> *next_it = new_it;
        	C_it<T> *last_it = new_it+n;
		new_it->itr_->prev=it.itr_->prev;
		while (next_it!=last_it){
			new_it->itr_->value=val;
			new_it++;
			next_it->itr_->next=new_it->itr_;
			new_it->itr_->prev=next_it->itr_;
			next_it++;
		}
	       	new_it->itr_->next=it.itr_;
		it.itr_->prev->next=new_it->itr_;
		it.itr_->prev=new_it->itr_;
		return *first_it;
	}

	iterator begin(void) const
	{
		return origin_;
	}

	const_iterator cbegin(void) const
	{
		return origin_;
	}
};

#endif
