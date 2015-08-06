#ifndef _SIMPLE_TYPES_H_
#define _SIMPLE_TYPES_H_

#include "key.h"
#include "basic_types.h"


class GENPROB :  public REAL {
public:
	GENPROB();
};

class FREQ :  public REAL {
public:
	FREQ();
};

class ERROR :  public REAL {
public:
	ERROR();
};

class LOGLIKE :  public REAL {
public:
	LOGLIKE();
};

class SMPNAME :  public SIZE {
public:
	SMPNAME();
};

class ROWID :  public SIZE {
public:
	ROWID();
};

class COV :  public COUNT {
public:
	COV();
};

class SMPSIZE :  public SIZE {
public:
	SMPSIZE();
};

class SOMATIC :  public BOOL {
public:
	SOMATIC();
};

class VALID :  public BOOL {
public:
	VALID();
};

class STRBIAS :  public REAL {
public:
	STRBIAS();
};

class SB :  public REAL {
public:
	SB();
};

class MAPQ0 :  public COUNT {
public:
	MAPQ0();
};

class HAPMAP2 :  public BOOL {
public:
	HAPMAP2();
};

class END :  public SIZE {
public:
	END();
};

class REF :  public GENOTYP {
public:
	REF();
};

class ANCTYPE :  public GENOTYP {
public:
	ANCTYPE();
};

class NCUT :  public SIZE {
public:
	NCUT();
};

class GOF :  public REAL {
public:
	GOF();
};

class DEFAULT :  public BOOL {
public:
	DEFAULT();
};

class POLYLL :  public REAL {
public:
	POLYLL();
};

class HWELL :  public REAL {
public:
	HWELL();
};

class MAXLL :  public REAL {
public:
	MAXLL();
};

class EFCHROM :  public REAL {
public:
	EFCHROM();
};
#endif
