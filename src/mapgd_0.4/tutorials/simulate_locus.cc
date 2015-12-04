#include "../datatypes/locus.h"
#include "../typedef.h"
 
int main(int argc, char* argv[])	// 
{
	uint64_t sum_cov=0;			// an unsigned 64 bit integer
	locus l(20);
	for (int x=0; x<10; ++x){
        std::cout << l << '\n';	
	}
	return 0;
}
