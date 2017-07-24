#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "typedef.h"

template <class R, class T>
class Constants{
private:
	size_t _size;
	void (**cfn) (Constants <R, T> &, const T &);
public:
	R* c;
	Constants(const size_t &size, void (**fn) (Constants <R, T> &, const T &) )
	{
		_size=size;
		c=(R *) calloc( size, sizeof(R) );
		cfn=(void (**) (Constants <R, T> &,const T &)) calloc( size, sizeof( cfn ));
		memcpy(cfn, fn, sizeof( cfn )*_size);
	}
	~Constants(void)
	{
		free(c);
		free(cfn);
	}
	void recalculate(const T &t)
	{
		for (int x=0; x<_size; x++) (*cfn[x])(*this, t);
	}
};
#endif
