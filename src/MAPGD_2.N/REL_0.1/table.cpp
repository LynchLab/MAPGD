//currently we are using a tab seperated file to store data, this should be changed to a binary file soon.

#include  "readtable.h"

using namespace std;

ll_t max(ll_t a, ll_t b){
	if (a>b) return a;
	else return b;
}; 

void streamtable(void){
	streamtable("", "");
};

void streamtable(const char *inname, const char *outname){
	istream *in;
	ostream *out;

	ifstream infile=NULL;
	ifstream outfile=NULL;

	if (inname!="") {
		infile.open(inname);
		if (!in.is_open()) {cerr << "cannot open " << inname << endl; exit(1);};
	}
	if (outname!="") {
		outfile.open(outname, ios::binary);
		if (!out.is_open()) {cerr << "cannot open " << outname << endl; exit(1);};
	{

	//our first order of buisness is to figure out the number of coloumns in the file;
	float_t MM, Mm, mm, P;
	size_t N;

	size_t SIZE;
	bool first=true;

	POPGL pgl;

	std::vector <std::string> args;
	std::vector <std::string>::iterator it;
	std::vector <std::string>::iterator end;

	readheader();

	while (getline(in, buffer)!=EOF){
		args=split(line, delim_column);
		end=ars.end();
		if (args.size()>4) 
		P=atof(ars[?].c_str());

		pgl.P=P;
		while (it!=end ){
			MM=atof(it->c_str()); ++t; if (it==end ) { cerr << "error parsing " << inname << " : " << buffer << endl; exit(1);};
			Mm=atof(it->c_str()); ++t; if (it==end ) { cerr << "error parsing " << inname << " : " << buffer << endl; exit(1);};
			mm=atof(it->c_str()); ++t; if (it==end ) { cerr << "error parsing " << inname << " : " << buffer << endl; exit(1);};
			N=atoi(it->c_str()); ++t;
			pgl.add(MM, Mm, mm, N);
		}
		if (first) {
			SIZE=pgl.gl.size();
			out.write((char *)&SIZE,sizeof(size_t) ); 
			first=false;
		}
		*out.write((char *)&pgl.P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; x++){
			out.write((char *)&pgl.gl[x].lMM,sizeof(ll_t) );
			out.write((char *)&pgl.gl[x].lMm,sizeof(ll_t) );
			out.write((char *)&pgl.gl[x].lmm,sizeof(ll_t) );
			out.write((char *)&pgl.gl[x].N,sizeof(size_t) );
		};
		pgl.clear(); 
	};
}


table_t readtable(const char *filename){
	ifstream file;
	string buffer;
	string token;
	file.open(filename);
	table_t table;
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};

	//our first order of buisness is to figure out the number of coloumns in the file;
	ll_t MM, Mm, mm, P;
	size_t N;
	POPGL pgl;
	while (!file.eof()){
		getline(file,buffer);
		if (file.eof() ) break;
   		boost::tokenizer<> tok(buffer);
		boost::char_separator<char> sep("\t");
		boost::tokenizer<boost::char_separator<char> > tokens(buffer, sep);
		boost::tokenizer<boost::char_separator<char> >::iterator t=tokens.begin();
		for (int x=0; x<4; ++x) ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
		P=atof(t->c_str());
		pgl.P=P;
		for (int x=0; x<2; ++x) ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
		while (t!=tokens.end() ){
			MM=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			Mm=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			mm=atof(t->c_str()); ++t; if (t==tokens.end() ) { cerr << "error parsing " << filename << endl; exit(1);};
			N=atoi(t->c_str()); ++t;
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(POPGL(pgl));
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

	ll_t MM, Mm, mm;
	size_t N, SIZE;
	file.read((char *)&SIZE,sizeof(size_t) );
	POPGL pgl;
	while (!file.eof()){
		file.read((char *)&pgl.P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(ll_t) );
			file.read((char *)&Mm,sizeof(ll_t) );
			file.read((char *)&mm,sizeof(ll_t) );
			file.read((char *)&N,sizeof(size_t) );
			pgl.add(MM, Mm, mm, N);
		}
		table.push_back(POPGL(pgl));
		pgl.clear(); 
	};
	return table;
};

std::map<PAIRGL, size_t> readcounts(const char *filename, size_t a, size_t b, size_t s){
	std::map<PAIRGL, size_t> counts;
	ifstream file;
	file.open(filename, ios::binary);
	
	if (!file.is_open()) {cerr << "cannot open " << filename << endl; exit(1);};


	ll_t MM, Mm, mm, AN=0, BN=0, tN=0, sampled=0;
	size_t SIZE, N;
	file.read((char *)&SIZE,sizeof(size_t) );
	POPGL pgl;
	while (!file.eof()){
		file.read((char *)&pgl.P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.read((char *)&MM,sizeof(ll_t) );
			file.read((char *)&Mm,sizeof(ll_t) );
			file.read((char *)&mm,sizeof(ll_t) );
			file.read((char *)&N,sizeof(size_t) );
			if (x==a) AN+=N;
			if (x==b) BN+=N;
		}
		tN+=1;
	};
	file.close();
	file.open(filename, ios::binary);

        file.read((char *)&SIZE,sizeof(size_t) );
	ll_t amin=max(AN/tN/2., 1), amax=AN*2./tN, bmin=max(BN/tN/2., 1), bmax=BN*2./tN;
	if (amax<6 || bmax<6){
		cout << "Sampled: " << sampled << "(" << a << ":" << amin << "-" << amax << ", " << b << ":" << bmin << "-" << bmax << ") ";
		return counts;
	}
        while (!file.eof()){
                file.read((char *)&pgl.P,sizeof(ll_t) );
                for (size_t x=0; x<SIZE; ++x){
                        file.read((char *)&MM,sizeof(ll_t) );
                        file.read((char *)&Mm,sizeof(ll_t) );
                        file.read((char *)&mm,sizeof(ll_t) );
                        file.read((char *)&N,sizeof(size_t) );
                        pgl.add(MM, Mm, mm, N);
                }
                if (pgl.gl[a].N>amin && pgl.gl[a].N<amax && pgl.gl[b].N>bmin && pgl.gl[b].N<bmax){sampled++; counts[PAIRGL (pgl.gl[a].lMM, pgl.gl[a].lMm, pgl.gl[a].lmm, pgl.gl[b].lMM, pgl.gl[b].lMm, pgl.gl[b].lmm, pgl.P)]+=1;}
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
	POPGL pgl;
	while (it!=end){
		file.write((char *)&it->P,sizeof(ll_t) );
		for (size_t x=0; x<SIZE; ++x){
			file.write((char *)&it->gl[x].lMM,sizeof(ll_t) );
			file.write((char *)&it->gl[x].lMm,sizeof(ll_t) );
			file.write((char *)&it->gl[x].lmm,sizeof(ll_t) );
			file.write((char *)&it->gl[x].N,sizeof(size_t) );
		}
		it++;
	};
};
