#ifndef PROFILE_HPP_
#define PROFILE_HPP_	1

#include <cstdio>
#include <cstring>
#include <iostream>
#include <climits>

class  profile{
	private:
		bool am_open;
	public:
	static char empty[30];
	profile();
	FILE *instream;
	char id1[30];
	char id2[30];
	int n[5];
	profile* open(const char*, const char *);
	bool is_open(void);
	int read(void);
	void close(void);
};

void sync(int, profile*);

#endif
