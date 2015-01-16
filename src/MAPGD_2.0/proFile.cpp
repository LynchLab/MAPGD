#include "proFile.h"

/** @brief Read a .pro file in the five column format.
  * @returns 0 iff successful
**/

const std::string profile::names_ ="ACGT*";

int profile::read(void){
	return profile::read(NONE);
};

int profile::read(int arg){
	if (donothing_) {donothing_=false; return 0;};
	std::string line;
	sorted_[0]=0;
	sorted_[1]=1;
	sorted_[2]=2;
	sorted_[3]=3;
	sorted_[4]=4;
	std::vector <std::string> fields;
	if (std::getline(*in, line)!=NULL){
		//READ A PRO FILE
		if (columns_==5){
			/* check for start of new scaffold*/
			if (line[0]!='>'){
				fields=split(line, delim);
				site_.id1=fields[0];
				if (fields.size()-1/4!=samples_){
					if (masked_!=NULL) delete masked_; 
					samples_=(fields.size()-1)/4;
					masked_=new bool[samples_];
				};
				while (site_.sample.size()<samples_) site_.sample.push_back(quartet_t() );
				for (unsigned int x=0; x<samples_; ++x){
					site_.sample[x].base[0]=atoi(fields[1+x*4].c_str() );
					site_.sample[x].base[1]=atoi(fields[2+x*4].c_str() );
					site_.sample[x].base[2]=atoi(fields[3+x*4].c_str() );
					site_.sample[x].base[3]=atoi(fields[4+x*4].c_str() );
				};
			} else {
				line.erase(line.begin());
				site_.id0=line;
				read(arg);
			}
		} else if (columns_==6){
			fields=split(line, delim);
			site_.id0=fields[0];
			site_.id1=fields[1];
			samples_=fields.size()-2;
			while (site_.sample.size()<samples_) site_.sample.push_back(quartet_t());
			for (unsigned int x=0; x<samples_; ++x){
				site_.sample[x].base[0]=atoi(fields[2+x*4].c_str() );
				site_.sample[x].base[1]=atoi(fields[3+x*4].c_str() );
				site_.sample[x].base[2]=atoi(fields[4+x*4].c_str() );
				site_.sample[x].base[3]=atoi(fields[5+x*4].c_str() );
			};
		}
		if (arg==PEAK) donothing_=true;
		return 0;
	}
	return EOF;
};

/** @brief opens a .pro file in the modes "r" or "w".
  * @returns a pointer to the profile
**/
profile* profile::open(const char* filename, const char mode){
	sorted_[0]=0;
	sorted_[1]=1;
	sorted_[2]=2;
	sorted_[3]=3;
	sorted_[4]=4;
	open_=false;

	switch ( mode ) {
		case 'r':
			inFile.open(filename, std::ifstream::in);
			if (!inFile){
				std::cerr << "cannot open " << filename << " for reading." << std::endl;				
				exit(0);
			};
			in=&inFile;
			break;
		case 'w':
			outFile.open(filename, std::ofstream::out);
			if (!outFile){
				std::cerr << "cannot open " << filename << " for reading." << std::endl;				
				exit(0);
			};
			out=&outFile;
			break;
		default :
			std::cerr << "unkown filemode " << std::endl;
			exit(0);
	
	}
	columns_=5;
	samples_=0;
	open_=true;
	masked_=NULL;
	donothing_=false;
	return this;
}

profile* profile::open(char mode)
{
	sorted_[0]=0;
	sorted_[1]=1;
	sorted_[2]=2;
	sorted_[3]=3;
	sorted_[4]=4;
	open_=false;
	
	switch ( mode ) {
		case 'r':
			in=&std::cin;
			break;
		case 'w':
			out=&std::cout;
			break;
		default :
			std::cerr << "unkown filemode " << std::endl;
			exit(0);
	
	};
	columns_=5;
	samples_=0;
	masked_=NULL;
	open_=true;
	donothing_=false;
	return this;
}
/** @brief closes a .pro file and unsets members.
  * @no return value
**/
void profile::close(void){
	open_=false;
};

profile::profile(){
	open_=false;
};

bool profile::is_open(void){
	return open_;
}

void profile::sort(void){
	count_t total[5]={0};
	for (unsigned int s=0; s<samples_;++s){
		if (masked_[s]) continue;
		total[0]+=site_.sample[s].base[0];
		total[1]+=site_.sample[s].base[1];
		total[2]+=site_.sample[s].base[2];
		total[3]+=site_.sample[s].base[3];
		total[4]+=site_.sample[s].base[4];
	};
	if (total[sorted_[0]]<total[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (total[sorted_[1]]<total[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (total[sorted_[2]]<total[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (total[sorted_[0]]<total[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (total[sorted_[1]]<total[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
};

void profile::sort(count_t s){
	if (site_.sample[s].base[sorted_[0]]<site_.sample[s].base[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (site_.sample[s].base[sorted_[1]]<site_.sample[s].base[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (site_.sample[s].base[sorted_[2]]<site_.sample[s].base[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (site_.sample[s].base[sorted_[0]]<site_.sample[s].base[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (site_.sample[s].base[sorted_[1]]<site_.sample[s].base[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
};

count_t profile::getcount(count_t s, count_t c)
{
	if (c<5) return site_.sample[s].base[sorted_[c]];
	else return 0;
};

count_t profile::getcount(count_t c)
{
	count_t total=0;
	if (c<5){
		for (unsigned int s=0; s<samples_;++s){
			if (masked_[s]) continue;
			total+=site_.sample[s].base[sorted_[c]];
		}
		return total;
	}
	else return 0;
};

count_t profile::getcoverage(count_t s)
{
	if (s<samples_)	return site_.sample[s].base[0]+
				site_.sample[s].base[1]+
				site_.sample[s].base[2]+
				site_.sample[s].base[3];
	else {
		std::cerr << "Attempted to access a population that doesn't exist. Exiting." << std::endl;
		exit(0);
	};
};

count_t profile::getcoverage()
{
	count_t total=0;
	for (unsigned int s=0; s<samples_;++s){
		if (masked_[s]) continue;
		total+=site_.sample[s].base[0]+
			site_.sample[s].base[1]+
			site_.sample[s].base[2]+
			site_.sample[s].base[3]+
			site_.sample[s].base[4];
	};
	return total;
};

name_t profile::getname(count_t c)
{
	if (c<5) return  profile::names_[sorted_[c]];
	else return '*';
};

std::string profile::getids(void)
{
	std::string str=site_.id0+'\t'+site_.id1;
	return str;
};

count_t profile::size(void){
	return samples_;
};

void profile::maskall(void){
	for (unsigned int s=0; s<samples_;++s){
		masked_[s]=true;
	};
};

void profile::unmask(count_t a){
	masked_[a]=false;
};
