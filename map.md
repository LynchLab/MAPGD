An example of how to use mapgd to get the average coverage of a genome. In this tutorial we\'ll see how to open a .pro file,
which is mapgd\'s way of representing sequence information, as well as how to do a little something with the information in 
.pro files. For more information on .pro files, please reffer to the developer's documation of mapgd.


```c++
#include "profile.h"					// For the profile class.
#include "Locus.h"					// For the Locus class.
#include "stdio.h"					// std::cout. If you are new to see, the 'std::cout' notation specifies that we are 
							// going to use member 'cout' of of the namespace 'std'. Namespaces help prevent an 
							// ? in complex programs. If you don't want to type std::this and std::that all of the 
							// time you can type 'using namespace std;'. However, I find it easier ot just type std::this, 
							// since it is six keystokes. 


int main (int argc, char *argv [])
{
	profile pro;					// Declare the profile.
	if (argc!=2){					// Check to see if the user passed the right number of arguments.
		std::cerr << "Usage foo bar.\n";	// If not print a totally useless error message,
		exit(0);				// and exit. (Hint: your error message shouldn't be totally useless,
							// only tutorials are allowed to do that.)
	}

	pro.open(argv[1], std::fstream::in);		// Open the profile in input mode.	

	if ( not (pro.is_open() ) ) {			// Check to see if the .pro file opened successfully. 
		std::cerr << "Oops!\n";			// If not print a totally useless error message,
		exit(0);				// and exit. 
	}

	Locus buffer;					// Declare a locus as a buffer for lines from the .pro file.

	count_t sum					// count_t is the typedef of depth of coverage type infomation in mapgd. count_t,
							// rather than uint16_t or unsigned int, or int, or whatever, should be used so that
							// commands can be compiled on different computers, since int ...

	while (pro.getline(buffer) ){			// Ye good old while loop. The loop will exit when the .pro file reaches EOF.
		sum+=buffer.count();			// Sum up the coverage from each site.
	};
	std::cout << double(sum)/double(pro.size())	// Print out the average coverage. 
	exit(0);					// FIN.
}
```
