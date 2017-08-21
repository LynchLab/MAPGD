#include "state.h"

const std::string State::file_name=".bgt";
const std::string State::table_name="STATE";
const bool State::binary=true;

const Registration State::registered=Registration(State::table_name, State::create);

State::State ()
{
	lz4_start_ = new char [LZ4_BUFFER_SIZE];
	lz4_ptr_=lz4_start_;
	lz4_last_ = lz4_start_;
	lz4_end_ = lz4_start_ + LZ4_BUFFER_SIZE;

	sites_=0;
	size_=0;
	cached_sites_=false;
}

void
State::set_stream(State_stream &stream) const
{
	stream.sites_=sites_;
	stream.size_=size_;
	if (cached_) stream.lz4_ptr=lz4_start_;
}

/*
State::State (const std::vector <std::string> &column)
{
	_sites=0;
	_size=0;
	compressed=false;
}*/
const std::string State::get_file_name(void) const {return this->file_name;}
const std::string State::get_table_name(void) const {return this->table_name;}

State::State(const std::vector <std::string> &column_names)
{
	//std::cerr << __FILE__ << ", " << __LINE__ << "const std::vector <std::string> &column_names called" << std::endl;
	if (column_names.size()==3) 
	{
		uint32_t sample_size=atoi(split(column_names[0], ':')[1].c_str() );
		uint32_t sites_size=atoi(split(column_names[1], ':')[1].c_str() );
		uint32_t buffer_size=atoi(split(column_names[2], ':')[1].c_str() );

		std::cerr << sample_size << ", " << sites_size << ", " << buffer_size << std::endl;
	
		size_=sample_size;
		sites_=sites_size;
		cached_sites_=sites_size;
		cached_=true;

		lz4_buffer_size_=LZ4_BUFFER_SIZE;

		lz4_start_ = new char [lz4_buffer_size_];
		lz4_ptr_=lz4_start_;
		lz4_last_ = lz4_start_+buffer_size;
		lz4_end_ = lz4_start_ + lz4_buffer_size_;
		block_size_=sizeof(uint32_t)*size_;
	}
}

State::State(const uint32_t &sample_size, const uint32_t &sites_size, const uint32_t &buffer_size)
{
	//std::cerr << __FILE__ << ", " << __LINE__ << "t uint32_t &sample_size, const uint32_t &sites_size, const uint32_t &buffer_size" << std::endl;
	lz4_buffer_size_=LZ4_BUFFER_SIZE;
	lz4_start_ = new char [lz4_buffer_size_];
	lz4_ptr_=lz4_start_;
	lz4_last_ = lz4_start_ + buffer_size;
	lz4_end_ = lz4_start_ + lz4_buffer_size_;

	size_=sample_size;
	sites_=sites_size;
	cached_sites_=sites_size;
	block_size_=sizeof(uint32_t)*size_;
	cached_=true;
}

State::State(const uint32_t &sample_size)
{
	//std::cerr << __FILE__ << ", " << __LINE__ << "const uint32_t &sample_size" << std::endl;
	size_=sample_size;
	sites_=0;
	cached_=false;
	cached_sites_=0;
	lz4_buffer_size_=LZ4_BUFFER_SIZE;

	lz4_start_ = new char [lz4_buffer_size_];
	lz4_ptr_=lz4_start_;
	lz4_last_ = lz4_start_;
	lz4_end_ = lz4_start_ + lz4_buffer_size_;
	block_size_=sizeof(uint32_t)*size_;
}

State::~State()
{
	if (lz4_start_){
		clear();
	//	std::cerr << "deleting " << size_t(lz4_start_) << std::endl;
		delete [] lz4_start_;
	}
}


void 
State::uncompress (uint32_t *a, uint32_t *b)
{
	if (!cached_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from uncached stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	if (sites_-- ) 
	{
		int ret;
		//std::cerr << (long int)(lz4_end_-lz4_ptr_) << ", " << (long int)(lz4_ptr_-lz4_start_) << ", " << (long int)(lz4_last_-lz4_ptr_) << ", " << (long int)(lz4_end_) << std::endl;
#ifndef NOLZ4
		ret=LZ4_decompress_fast (lz4_ptr_, (char *)a, block_size_);
#else
#endif
		if (ret > 0)
		{
			lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
#ifndef NOLZ4
		ret=LZ4_decompress_fast (lz4_ptr_, (char *)b, block_size_);	
#else
#endif
		if (ret > 0)
		{
			lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::uncompress (uint32_t *a, uint32_t *b, State_stream &stream) const
{
	if (!cached_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from uncached stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	if (stream.sites_-- ) 
	{
		int ret;
		//std::cerr << (long int)(lz4_end_-lz4_ptr_) << ", " << (long int)(lz4_ptr_-lz4_start_) << ", " << (long int)(lz4_last_-lz4_ptr_) << ", " << (long int)(lz4_end_) << std::endl;
#ifndef NOLZ4
		ret=LZ4_decompress_fast (stream.lz4_ptr, (char *)a, block_size_);
#else
#endif
		if (ret > 0)
		{
			stream.lz4_ptr+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
#ifndef NOLZ4
		ret=LZ4_decompress_fast (stream.lz4_ptr, (char *)b, block_size_);	
#else
#endif
		if (ret > 0)
		{
			stream.lz4_ptr+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::uncompress (uint32_t *a, uint32_t *b, const uint32_t &k)
{
	while (cached_sites_-sites_-- < k) 
	{
#ifndef NOLZ4
		lz4_ptr_+=LZ4_decompress_fast (lz4_ptr_, (char *)a, block_size_);
		lz4_ptr_+=LZ4_decompress_fast (lz4_ptr_, (char *)b, block_size_);
#else
#endif
	} 
	if (sites_==0) {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::cache (void)
{
	if(!cached_)
	{
		lz4_last_=lz4_ptr_;
		lz4_ptr_=lz4_start_;
		cached_=true;
		cached_sites_=sites_;
	}
}

void
State::rewind(void)
{
	if (cached_)
	{
		lz4_ptr_=lz4_start_;
		cached_=true;
		sites_=cached_sites_;
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to rewind uncached stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::increase_buffer_(void)
{
	char *new_lz4=new char[lz4_buffer_size_+LZ4_BUFFER_SIZE];
	memcpy(new_lz4, lz4_start_, lz4_buffer_size_);
	lz4_buffer_size_+=LZ4_BUFFER_SIZE;
	lz4_ptr_=new_lz4+(lz4_ptr_-lz4_start_);
	lz4_last_=new_lz4+(lz4_last_-lz4_start_);
	lz4_end_=new_lz4+lz4_buffer_size_;
	delete [] lz4_start_;
	lz4_start_=new_lz4;
}

void 
State::compress (const uint32_t *a, const uint32_t *b)
{
	int size;
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) increase_buffer_();
#ifndef NOLZ4
	size=LZ4_compress_default( (const char*) a, lz4_ptr_, block_size_,  size_t (lz4_end_-lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-lz4_ptr_) );
#else
#endif
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	lz4_ptr_+=size;
#ifndef NOLZ4
	size=LZ4_compress_default( (const char*) b, lz4_ptr_, block_size_,  size_t (lz4_end_-lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-lz4_ptr_) );
#else
#endif
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	lz4_ptr_+=size;
	++sites_;
	cached_=false;
	if (lz4_end_==lz4_ptr_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}


void 
State::clear(void)
{
	lz4_ptr_=lz4_start_;
	lz4_last_=lz4_start_;
	cached_=false;
	sites_=0;
	cached_sites_=0;
}

static uint32_t mask[32]={0x00000001,	0x00000002,	0x00000004,	0x00000008,
					0x00000010,	0x00000020,	0x00000040,	0x00000080,
					0x00000100,	0x00000200,	0x00000400, 	0x00000800,	
					0x00001000, 	0x00002000, 	0x00004000, 	0x00008000,
					0x00010000,	0x00020000,	0x00040000,	0x00080000,
					0x00100000,	0x00200000,	0x00400000,	0x00800000,
					0x01000000,	0x02000000,	0x04000000, 	0x08000000,	
					0x10000000, 	0x20000000, 	0x40000000, 	0x80000000};
void 
State::write (std::ostream& out) const 
{
	uint32_t *set0=new uint32_t[size_];
	uint32_t *set1=new uint32_t[size_];

	State_stream stream;
	set_stream(stream);

	while (stream.sites_>0)
	{
		uncompress(set0, set1, stream);
		for (size_t b=0; b<32; ++b) 
		{
			for (size_t y=0; y<size_; ++y)
			{
				if (set0[y] & mask[b]) {
					if ( (set1[y] & mask[b] ) ) 
					{
						out << "	1	1";
					} else {
						out << "	1	0";
					}
				} else {
					if ( (set1[y] & mask[b]) )
					{
						out << "	0	1";
					} else {
						out << "	0	0";
					}
				}
			}
			out << std::endl;
		}
	}
	delete set0;
	delete set1;
}
	
double 
State::compression_ratio (void) const
{
	if (cached_)
	{
		return double(lz4_last_-lz4_start_)/double(2*block_size_*sites_);
	}
	else return 0;
}

void 
State::read (std::istream& in)
{
	uint32_t *set0=new uint32_t[size_];
	uint32_t *set1=new uint32_t[size_];

	uint32_t w, z, sites=sites_;
	std::cerr << "Sites:" << sites_ << std::endl;
	clear();
	for (size_t x=0; x<sites; ++x)
	{
		memset(set0, 0, sizeof(uint32_t)*size_);
		memset(set1, 0, sizeof(uint32_t)*size_);
		for (size_t b=0; b<32; ++b) 
		{
			for (size_t y=0; y<size_; ++y)
			{
				in >> w;
				in >> z;
				if (w) set0[y]+=mask[b];
				if (z) set1[y]+=mask[b];
			}
		}
		compress(set0, set1);
	}
	delete set0;
	delete set1;
	cache();
}

/*
void
State::transpose(void)
{
	uint32_t *set0_vert=new uint32_t[size_];
	uint32_t *set1_vert=new uint32_t[size_];
	
	char *set0_horz=new char[size_*4];
	char *set1_horz=new char[size_*4];

	row_size=size_ 
	padding=
	per_row=
	jM
	jR
	for (size_t j=0; j<size_; ++j)
	{
		for (size_t i; i<32; ++i)
		{
			set0_horz[jM+i*per_row] |= (set0_vert[j] & mask[i] != 0) << jR;
			set1_horz[jM+i*per_row] |= (set0_vert[j] & mask[i] != 0) << jR;
		}
	}
}*/

std::string 
State::header (void) const
{
	return "@NS:"+std::to_string(size_)+"\tGS:"+std::to_string(sites_)+"\tBS:"+std::to_string(lz4_last_-lz4_start_)+'\n';
	//return double(lz4_last_-lz4_start_)/double(2*block_size_*sites_);
}

void 
State::read_binary (std::istream& in)
{
	in.read(lz4_start_, lz4_last_-lz4_start_ );
	cached_=true;
}

void 
State::write_binary (std::ostream& out) const 
{
	out.write(lz4_start_, lz4_last_-lz4_start_ );
}
