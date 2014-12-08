#include "ReadFile.hpp"

int profile::read(void){
	//READ A PRO FILE
	char c = fgetc (instream);
	fseek(instream, -1, SEEK_CUR);

	/* check for start of new scaffold*/

	if (c!='>'){if (fscanf(instream,"%s\t%i\t%i\t%i\t%i", id2, &n[1], &n[2], &n[3], &n[4])==EOF) return EOF;}
	else {
		if (fscanf(instream, "%s", id1)==EOF) return EOF;
		read();
	}
	return 0;
};

void profile::open(const char* filename, const char* mode){
	instream = fopen(filename, mode);
	if(instream==NULL) {std::cerr << "failed to open file " << filename << " for reading\n"; exit(0); }
};

void sync(profile*);
