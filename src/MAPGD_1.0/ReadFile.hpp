#ifndef PROFILE_HPP_
#define PROFILE_HPP_	1

#include <cstdio>
#include <iostream>

class  profile{
	private:
	public:
	FILE *instream;
	char id1[30];
	char id2[30];
	int n[5];
	void open(const char*, const char *);
	int read(void);
};

void sync(profile*);

#endif
