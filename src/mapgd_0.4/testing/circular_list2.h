#ifndef _CICULAR_LIST_H_
#define _CICULAR_LIST_H_

#include <iostream>
#include <list>

template <class Type>
class Circular_list : public std::list<Type> {
public:
	using std::list<Type>::iterator;
	using std::list<Type>::const_iterator;
	Circular_list (const Type &T)
	{
		
	}
};

#endif
