#ifndef _BASE_CLASSES_H_
#define _BASE_CLASSES_H_

//#include "quartet.h"
//#include "locus.h"
//#include "allele.h"
#include "../typedef.h"

/// a class that handels the formating of data.
/* It also makes sure that any data being passed through map_pipe has all the 
 * commands neccisary to format the data properly for reading and writting.
 */
class data{
protected:	
	char keyname_[8];
	uint8_t keynum_;
	std::string keydesc_;
public:
	virtual size_t sizeb()=0;
	virtual size_t sizet()=0;

	virtual char * to_text()=0;
	virtual void from_text(char *)=0;

	const char * get_name(void) const;
	const uint8_t get_num(void) const;
	const std::string get_desc(void) const;
};


class CHAR7 : public data {
public:
	size_t sizeb(void){return 8;};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
};

class CHAR13 : public data {
public:
	size_t sizeb(void){return 14;};
	size_t sizet(void){return 13;};
	char * to_text();
	void from_text(char *);
};

class UINT8 : public data {
public:
	size_t sizeb(void){return 1;};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
};

/// every key has a name.
class KEYNAME : public CHAR7 {
public:
	KEYNAME();
	KEYNAME(const char *this_keyname);
};

/// every key has a number.
class KEYNUM : public UINT8 {
public:
	KEYNUM(void);
	KEYNUM(uint8_t keynum);
};

/// every key has a description.
class KEYDESC : public data {
public:
	size_t sizeb(void){return 0;};
	size_t sizet(void){return 0;};
	char * to_text() {return 0;};
	void from_text(char *) {};
	KEYDESC();
	KEYDESC(std::string desc);
};


class SIZE : public data {
public:
	size_t sizeb(void){return sizeof(size_t);};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
};

class REAL : public data{
public:
	size_t sizeb(void){return sizeof(long double);};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
};

class GENOTYP : public REAL {
public:
	size_t sizeb(void){return sizeof(gt_t);};
	size_t sizet(void){return 7;};
};

class COUNT : public SIZE {
public:
	size_t sizeb(void){return sizeof(count_t);};
	size_t sizet(void){return 7;};
};

class BOOL : public data {
public:
	size_t sizeb(void){return sizeof(bool);};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
};

class CHROM : public data {
public:
	size_t sizeb(void){return sizeof(char)*7;};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
	CHROM();
};

class POS : public data {
public:
	size_t sizeb(void){return sizeof(id1_t);};
	size_t sizet(void){return 7;};
	char * to_text();
	void from_text(char *);
	POS();
};
/*
class QUARTET :  public data {
public:
	size_t sizeb(void *){return sizeof(quartet_t);};
	size_t sizet(void *){return 14;};
	char * to_text();
	void from_text(char *);
	QUARTET();
};

class LOCUS :  public data {
public:
	size_t sizeb (void *source){ return 0;};//( (locus *) source )->size(); };
	size_t sizet (void *source){ return 0;};//( (locus *) source )->size()*8; };
	char * to_text();
	void from_text(char *);
	LOCUS();
};

class ALLELE :  public data {
public:
	size_t sizeb(void *){return sizeof(allele_t);};
	size_t sizet(void *){return 40;};
	char * to_text();
	void from_text(char *);
	ALLELE();
};*/
#endif
