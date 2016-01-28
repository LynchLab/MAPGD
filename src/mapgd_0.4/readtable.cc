//currently we are using a tab seperated file to store data, this should be changed to a binary file soon.

#include <fstream>
#include <cstdio>
#include <iostream>

#include  "readtable.h"

using namespace std;

const unsigned char MIN_QUAL = '!';
float_t max(float_t a, float_t b){
	if (a>b) return a;
	else return b;
}; 

void streamtable(const char *inname, const char *outname){
	ifstream in;
	ofstream out;

	in.open(inname);
	out.open(outname, ios::binary);

	if (!in.is_open()) {cerr << "cannot open " << inname << endl; exit(1);};
//		for (?) { cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << inname << " : " << buffer << endl; exit(1);};
	if (!out.is_open()) {cerr << "cannot open " << outname << endl; exit(1);};
//		for (?) { cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << inname << " : " << buffer << endl; exit(1);};

	//our first order of buisness is to figure out the number of coloumns in the file;
	float_t MM, Mm, mm, m;
	size_t N;

	size_t SIZE;
	bool first=true;

        std::string buffer;
        std::vector <std::string> column;
	std::vector <std::string> G;
	// args, arg;

	population_genotypes pgl;
	getline(in, buffer);
	
	column=split(buffer, '\t');
	std::vector <std::string> names(column.begin()+6, column.end() );
	SIZE=names.size();
	
	while (!in.eof()){
		getline(in, buffer);
		if (buffer[0]=='@') {
			continue;
		}
		if (in.eof() ) break;
		m=atof(column[5].c_str() ); //or 4!
		pgl.m=m;
		column=split(buffer, '\t');
		if (column.size()!=6+SIZE) {
			std::cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << inname << " : " << buffer << std::endl; 
			exit(0);
		}
		for (int x=6; x<6+SIZE; ++x){
			G=split(column[x], '/');
			if (G.size()==4) {
				MM=atof(G[0].c_str());
				Mm=atof(G[1].c_str());
				mm=atof(G[2].c_str());
				N=atoi(G[3].c_str());
				pgl.add(MM, Mm, mm, N);
			} else {
				std::cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << inname << " : " << buffer << std::endl; 
				exit(0);
			}
		}
		if (first) {
			SIZE=pgl.likelihoods.size();
			out.write((char *)&SIZE,sizeof(size_t) ); 
			first=false;
		}
		out.write((char *)&pgl.m,sizeof(float_t) );
		for (size_t x=0; x<SIZE; x++){
			out.write((char *)&pgl.likelihoods[x].lMM,sizeof(float_t) );
			out.write((char *)&pgl.likelihoods[x].lMm,sizeof(float_t) );
			out.write((char *)&pgl.likelihoods[x].lmm,sizeof(float_t) );
			out.write((char *)&pgl.likelihoods[x].N,sizeof(size_t) );
		};
		pgl.clear(); 
		buffer="";
	};
}


table_t readtable(const char *filename){
	ifstream file;

	file.open(filename);
	table_t table;
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	//our first order of buisness is to figure out the number of coloumns in the file;
	float_t MM, Mm, mm, m;

        std::string buffer;
        std::vector <std::string> column;
	std::vector <std::string> G;
	//args, arg;

	size_t N;
	population_genotypes pgl;

	column=split(buffer, '\t');
	std::vector <std::string> names(column.begin()+6, column.end() );
	int SIZE=names.size();

	while (!file.eof()){
		getline(file,buffer);
		if (file.eof() ) break;

		m=atof(column[5].c_str() );	//or 4!
		pgl.m=m;
		column=split(buffer, '\t');

		if (column.size()!=6+SIZE) {
			std::cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << filename << " : " << buffer << std::endl; 
			exit(0);
		}

		for (int x=6; x<6+SIZE; ++x){
			G=split(column[x], '/');
			if (G.size()==4) {
				MM=atof(G[0].c_str() );
				Mm=atof(G[1].c_str() );
				mm=atof(G[2].c_str() );
				N=atoi(G[3].c_str() );
				pgl.add(MM, Mm, mm, N);
			} else {
				std::cerr << __FILE__ << ":" << __LINE__ << ":" << " error parsing " << filename << " : " << buffer << std::endl; 
				exit(0);
			}
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(population_genotypes(pgl));
		pgl.clear(); 
		buffer="";
	};
	return table;
};

table_t readbin(const char *filename){
	ifstream file;
	string buffer;
	file.open(filename, ios::binary);
	table_t table;
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	float_t MM, Mm, mm;
	size_t N, SIZE;
	file.read((char *)&SIZE,sizeof(size_t) );
	population_genotypes pgl;
	while (!file.eof()){
		file.read((char *)&pgl.m,sizeof(float_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(float_t) );
			file.read((char *)&Mm,sizeof(float_t) );
			file.read((char *)&mm,sizeof(float_t) );
			file.read((char *)&N,sizeof(size_t) );
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(population_genotypes(pgl));
		pgl.clear(); 
	};
	return table;
};

std::map<genotype_pair, size_t> readcounts(const char *filename, size_t a, size_t b, size_t s){
	std::map<genotype_pair, size_t> counts;
	ifstream file;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};


	float_t MM, Mm, mm, AN=0, BN=0, tN=0, sampled=0;
	size_t SIZE, N;
	file.read((char *)&SIZE,sizeof(size_t) );
	population_genotypes pgl;
	while (!file.eof()){
		file.read((char *)&pgl.m,sizeof(float_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(float_t) );
			file.read((char *)&Mm,sizeof(float_t) );
			file.read((char *)&mm,sizeof(float_t) );
			file.read((char *)&N,sizeof(size_t) );
			if (x==a) AN+=N;
			if (x==b) BN+=N;
		}
		tN+=1;
	};
	file.close();
	file.open(filename, ios::binary);

        file.read((char *)&SIZE,sizeof(size_t) );
	float_t amin=max(AN/tN/2., 1), amax=AN*2./tN, bmin=max(BN/tN/2., 1), bmax=BN*2./tN;
	if (amax<6 || bmax<6){
		cout << "Sampled: " << sampled << "(" << a << ":" << amin << "-" << amax << ", " << b << ":" << bmin << "-" << bmax << ") ";
		return counts;
	}
        while (!file.eof()){
                file.read((char *)&pgl.m,sizeof(float_t) );
                for (size_t x=0; x<SIZE; ++x){
                        file.read((char *)&MM,sizeof(float_t) );
                        file.read((char *)&Mm,sizeof(float_t) );
                        file.read((char *)&mm,sizeof(float_t) );
                        file.read((char *)&N,sizeof(size_t) );
                        pgl.add(MM, Mm, mm, N);
                }
                if (pgl.likelihoods[a].N>amin && pgl.likelihoods[a].N<amax && pgl.likelihoods[b].N>bmin && pgl.likelihoods[b].N<bmax){sampled++; counts[genotype_pair (pgl.likelihoods[a].lMM, pgl.likelihoods[a].lMm, pgl.likelihoods[a].lmm, pgl.likelihoods[b].lMM, pgl.likelihoods[b].lMm, pgl.likelihoods[b].lmm, pgl.m)]+=1;}
                pgl.clear(); 
   		if (sampled==s && s!=0) break;
        };
	cout << "Sampled: " << sampled << ", (" << a << ":" << amin << "-" << amax << ", " << b << ":" << bmin << "-" << bmax << "), ";
	file.close();
	return counts;
};

void writebin(table_t table, const char *filename){
	ofstream file;
	string buffer;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	size_t SIZE;
	SIZE=table[1].size();

	table_t::iterator it=table.begin();
	table_t::iterator end=table.end();
	file.write((char *)&SIZE,sizeof(size_t) );
	population_genotypes pgl;
	while (it!=end){
		file.write((char *)&it->m,sizeof(float_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.write((char *)&it->likelihoods[x].lMM,sizeof(float_t) );
			file.write((char *)&it->likelihoods[x].lMm,sizeof(float_t) );
			file.write((char *)&it->likelihoods[x].lmm,sizeof(float_t) );
			file.write((char *)&it->likelihoods[x].N,sizeof(count_t) );
		}
		it++;
	};
};
