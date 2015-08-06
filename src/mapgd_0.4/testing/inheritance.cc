#include "../datatypes/key.h"
#include "../datatypes/basic_types.h"
#include "../typedef.h"

class FLAG : public BOOL {
public:
	const std::type_info& type(void ){return typeid(*this);};
	FLAG(){};
};

class GENPROB :  public REAL {
private:
        real_t value_;
public:
        const std::type_info& type(void ){return typeid(*this);};
        GENPROB(){};
};


int main () {
	FLAG flag;
	GENPROB genprob;
	return 0;
}
