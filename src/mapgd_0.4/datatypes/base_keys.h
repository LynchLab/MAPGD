#ifndef _BASE_KEYS_H_
#define _BASE_KEYS_H_

#include "key.h"
#include "formats.h"

class real_key : public key{
private:
	real_t *value_;
public:
	real_t *value();
};

#endif
