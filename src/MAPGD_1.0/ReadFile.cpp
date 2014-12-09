#include "ReadFile.hpp"

/** @brief Read a .pro file in the five column format.
  * @returns 0 iff successful
**/
int profile::read(void){
	//READ A PRO FILE
	char c = fgetc (instream);
	fseek(instream, -1, SEEK_CUR);

	/* check for start of new scaffold*/

	if (c!='>'){if (fscanf(instream,"%s\t%i\t%i\t%i\t%i", id2, &n[1], &n[2], &n[3], &n[4])==EOF){close(); return EOF;} }
	else {
		if (fscanf(instream, "%s", id1)==EOF){close(); return EOF;}
		read();
	}
	return 0;
};

/** @brief opens a .pro file in the modes "r" or "w".
  * @returns a pointer to the profile
**/
profile* profile::open(const char* filename, const char* mode){
	am_open=false;
	instream = fopen(filename, mode);
	if(instream==NULL) {std::cerr << "failed to open file " << filename << " for reading\n"; return NULL; }
	am_open=true;
};

/** @brief closes a .pro file and unsets members.
  * @no return value
**/
void profile::close(void){
	memset(id1, 0, sizeof(char)*30);
	memset(id2, 0, sizeof(char)*30);
	memset(n, 0, sizeof(int)*5);
	am_open=false;
};

profile::profile(){
	memset(id1, 0, sizeof(char)*30);
	memset(id2, 0, sizeof(char)*30);
//	memset(empty, 0, sizeof(char)*30);
	memset(n, 0, sizeof(int)*5);
	am_open=false;
};

bool profile::is_open(void){
	return am_open;
}
/** @brief Reads an array of .pro file one bp at a time.
  * @returns 0 iff successful
**/
void sync(int proc, profile *prof){
	//check to see that prof is initalized;
/*
	int minid1=INT_MAX;
	int minid2=INT_MAX;
	for (int x=0; x<proc; ++x){
		if (prof[x].is_open() ) if (strcmp(prof[x].id1, profile().empty)==0) if(prof[x].read()==EOF) prof[x].close(); 
		if (hash(prof[x].id1)<minid1) minid1=prof[x].id1;
		if (prof[x].id2<minid2) minid2=prof[x].id2;
	}
	for (int x=0; x<proc; ++x){
		if (prof[x].is_open() ){
			if (prof[x].id1==minid1 && prof[x].id2==minid2 && prof[x].read()==EOF) prof[x].close();
			else memset(prof[x].n, 0, sizeof(int)*5);
		};
	}
	return 0;*/
};
