#include "../datatypes/locus.h"
#include "../typedef.h"
 
int main(int argc, char* argv[])	// 
{
	uint64_t sum_cov=0;			// an unsigned 64 bit integer
	id1_t site_count=0;			// id1_t specifies	
	for (locus site(20);;) {
		std::cin >> site;
		if(!std::cin) break;
		sum_cov+=site.get_coverage();
		++site_count;
	}
	std::cout << sum_cov << "/" << site_count << std::endl;
	return 0;
}
