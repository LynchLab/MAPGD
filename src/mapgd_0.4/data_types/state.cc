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

	temp_lz4_ptr_=lz4_ptr_;

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
	transposed_=false;
	//std::cerr << __FILE__ << ", " << __LINE__ << "const std::vector <std::string> &column_names called" << std::endl;
	if (column_names.size()==3) 
	{
		uint32_t sample_size=atoi(split(column_names[0], ':')[1].c_str() );
		uint32_t sites_size=atoi(split(column_names[1], ':')[1].c_str() );
		uint32_t buffer_size=atoi(split(column_names[2], ':')[1].c_str() );

		//std::cerr << sample_size << ", " << sites_size << ", " << buffer_size << std::endl;
	
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

		temp_lz4_ptr_=lz4_ptr_;
	}
	k_=0;
}

State::State(const uint32_t &sample_size, const uint32_t &sites_size, const uint32_t &buffer_size)
{
	transposed_=false;
	//std::cerr << __FILE__ << ", " << __LINE__ << " t uint32_t &sample_size, const uint32_t &sites_size, const uint32_t &buffer_size" << std::endl;
	lz4_buffer_size_=LZ4_BUFFER_SIZE;
	lz4_start_ = new char [lz4_buffer_size_];
	lz4_ptr_=lz4_start_;
	lz4_last_ = lz4_start_ + buffer_size;
	lz4_end_ = lz4_start_ + lz4_buffer_size_;
	temp_lz4_ptr_=lz4_ptr_;

	size_=sample_size;
	sites_=sites_size;
	cached_sites_=sites_size;
	block_size_=sizeof(uint32_t)*size_;
	cached_=true;
	k_=0;
}

State::State(const uint32_t &sample_size)
{
	transposed_=false;
	//std::cerr << __FILE__ << ", " << __LINE__ << " const uint32_t &sample_size:" << size_t(this) << std::endl;
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
	temp_lz4_ptr_=lz4_ptr_;
	k_=0;
}

State::~State()
{
	if (lz4_start_){
		clear();
		//std::cerr << "deleting " << size_t(lz4_start_) << std::endl;
		delete [] lz4_start_;
		lz4_start_=NULL;
	}
}


void 
State::uncompress_inplace (uint32_t *a, uint32_t *b)
{
	temp_lz4_ptr_=lz4_ptr_;

	if (sites_) 
	{
		int ret;
		ret=LZ4_decompress_fast (temp_lz4_ptr_, (char *)a, block_size_);
		if (ret > 0)
		{
			temp_lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		ret=LZ4_decompress_fast (temp_lz4_ptr_, (char *)b, block_size_);	
		if (ret > 0)
		{
			temp_lz4_ptr_+=ret;
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
		ret=LZ4_decompress_fast (lz4_ptr_, (char *)a, block_size_);
		if (ret > 0)
		{
			lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		ret=LZ4_decompress_fast (lz4_ptr_, (char *)b, block_size_);	
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

//hackers_delight_transpose

/* Straight-line version of transpose32a & b. */

#define swap(a0, a1, j, m) t = (a0 ^ (a1 >> j)) & m; \
                           a0 = a0 ^ t; \
                           a1 = a1 ^ (t << j);

void transpose32c(const uint32_t A[32], uint32_t B[32]) {
   uint32_t m, t;
   uint32_t a0, a1, a2, a3, a4, a5, a6, a7,
            a8, a9, a10, a11, a12, a13, a14, a15,
            a16, a17, a18, a19, a20, a21, a22, a23,
            a24, a25, a26, a27, a28, a29, a30, a31;

   a0  = A[ 0];  a1  = A[ 1];  a2  = A[ 2];  a3  = A[ 3];
   a4  = A[ 4];  a5  = A[ 5];  a6  = A[ 6];  a7  = A[ 7];
   a8  = A[ 8];  a9  = A[ 9];  a10 = A[10];  a11 = A[11];
   a12 = A[12];  a13 = A[13];  a14 = A[14];  a15 = A[15];
   a16 = A[16];  a17 = A[17];  a18 = A[18];  a19 = A[19];
   a20 = A[20];  a21 = A[21];  a22 = A[22];  a23 = A[23];
   a24 = A[24];  a25 = A[25];  a26 = A[26];  a27 = A[27];
   a28 = A[28];  a29 = A[29];  a30 = A[30];  a31 = A[31];

   m = 0x0000FFFF;
   swap(a0,  a16, 16, m)
   swap(a1,  a17, 16, m)
   swap(a2,  a18, 16, m)
   swap(a3,  a19, 16, m)
   swap(a4,  a20, 16, m)
   swap(a5,  a21, 16, m)
   swap(a6,  a22, 16, m)
   swap(a7,  a23, 16, m)
   swap(a8,  a24, 16, m)
   swap(a9,  a25, 16, m)
   swap(a10, a26, 16, m)
   swap(a11, a27, 16, m)
   swap(a12, a28, 16, m)
   swap(a13, a29, 16, m)
   swap(a14, a30, 16, m)
   swap(a15, a31, 16, m)
   m = 0x00FF00FF;
   swap(a0,  a8,   8, m)
   swap(a1,  a9,   8, m)
   swap(a2,  a10,  8, m)
   swap(a3,  a11,  8, m)
   swap(a4,  a12,  8, m)
   swap(a5,  a13,  8, m)
   swap(a6,  a14,  8, m)
   swap(a7,  a15,  8, m)
   swap(a16, a24,  8, m)
   swap(a17, a25,  8, m)
   swap(a18, a26,  8, m)
   swap(a19, a27,  8, m)
   swap(a20, a28,  8, m)
   swap(a21, a29,  8, m)
   swap(a22, a30,  8, m)
   swap(a23, a31,  8, m)
   m = 0x0F0F0F0F;
   swap(a0,  a4,   4, m)
   swap(a1,  a5,   4, m)
   swap(a2,  a6,   4, m)
   swap(a3,  a7,   4, m)
   swap(a8,  a12,  4, m)
   swap(a9,  a13,  4, m)
   swap(a10, a14,  4, m)
   swap(a11, a15,  4, m)
   swap(a16, a20,  4, m)
   swap(a17, a21,  4, m)
   swap(a18, a22,  4, m)
   swap(a19, a23,  4, m)
   swap(a24, a28,  4, m)
   swap(a25, a29,  4, m)
   swap(a26, a30,  4, m)
   swap(a27, a31,  4, m)
   m = 0x33333333;
   swap(a0,  a2,   2, m)
   swap(a1,  a3,   2, m)
   swap(a4,  a6,   2, m)
   swap(a5,  a7,   2, m)
   swap(a8,  a10,  2, m)
   swap(a9,  a11,  2, m)
   swap(a12, a14,  2, m)
   swap(a13, a15,  2, m)
   swap(a16, a18,  2, m)
   swap(a17, a19,  2, m)
   swap(a20, a22,  2, m)
   swap(a21, a23,  2, m)
   swap(a24, a26,  2, m)
   swap(a25, a27,  2, m)
   swap(a28, a30,  2, m)
   swap(a29, a31,  2, m)
   m = 0x55555555;
   swap(a0,  a1,   1, m)
   swap(a2,  a3,   1, m)
   swap(a4,  a5,   1, m)
   swap(a6,  a7,   1, m)
   swap(a8,  a9,   1, m)
   swap(a10, a11,  1, m)
   swap(a12, a13,  1, m)
   swap(a14, a15,  1, m)
   swap(a16, a17,  1, m)
   swap(a18, a19,  1, m)
   swap(a20, a21,  1, m)
   swap(a22, a23,  1, m)
   swap(a24, a25,  1, m)
   swap(a26, a27,  1, m)
   swap(a28, a29,  1, m)
   swap(a30, a31,  1, m)

   B[ 0] = a0;   B[ 1] = a1;   B[ 2] = a2;   B[ 3] = a3;
   B[ 4] = a4;   B[ 5] = a5;   B[ 6] = a6;   B[ 7] = a7;
   B[ 8] = a8;   B[ 9] = a9;   B[10] = a10;  B[11] = a11;
   B[12] = a12;  B[13] = a13;  B[14] = a14;  B[15] = a15;
   B[16] = a16;  B[17] = a17;  B[18] = a18;  B[19] = a19;
   B[20] = a20;  B[21] = a21;  B[22] = a22;  B[23] = a23;
   B[24] = a24;  B[25] = a25;  B[26] = a26;  B[27] = a27;
   B[28] = a28;  B[29] = a29;  B[30] = a30;  B[31] = a31;
}

void
fast_trans(uint32_t const *A, uint32_t *B, int nrows, int ncols)
{
//	assert(nrows % 32 == 0 && ncols == 32)
	for (size_t x=0; x<nrows; x+=32)
		transpose32c(A+x, B+x);
}

//mischasan sse_transpose. 

#define II  i 

void
sse_trans(uint8_t const *inp, uint8_t *out, int nrows, int ncols)
{
#   define INP(x,y) inp[(x)*ncols/8 + (y)/8]
#   define OUT(x,y) out[(y)*nrows/8 + (x)/8]
    int rr, cc, i, h;
    union { __m128i x; uint8_t b[16]; } tmp;
//    assert(nrows % 8 == 0 && ncols % 8 == 0);

    // Do the main body in 16x8 blocks:
    for (rr = 0; rr <= nrows - 16; rr += 16) {
        for (cc = 0; cc < ncols; cc += 8) {
            for (i = 0; i < 16; ++i)
                tmp.b[i] = INP(rr + II, cc);
            for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
                *(uint16_t*)&OUT(rr,cc+II)= _mm_movemask_epi8(tmp.x);
        }
    }
    if (rr == nrows) return;

    // The remainder is a block of 8x(16n+8) bits (n may be 0).
    //  Do a PAIR of 8x8 blocks in each step:
    for (cc = 0; cc <= ncols - 16; cc += 16) {
        for (i = 0; i < 8; ++i) {
            tmp.b[i] = h = *(uint16_t const*)&INP(rr + II, cc);
            tmp.b[i + 8] = h >> 8;
        }
        for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1)) {
            OUT(rr, cc + II) = h = _mm_movemask_epi8(tmp.x);
            OUT(rr, cc + II + 8) = h >> 8;
        }
    }
    if (cc == ncols) return;

    //  Do the remaining 8x8 block:
    for (i = 0; i < 8; ++i)
        tmp.b[i] = INP(rr + II, cc);
    for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
        OUT(rr, cc + II) = _mm_movemask_epi8(tmp.x);
}

void
State::transpose()
{
	transposed_=!transposed_;
	uint32_t *a=new uint32_t [size_];
	uint32_t *b=new uint32_t [size_];

	uint32_t *c=new uint32_t [size_];
	uint32_t *d=new uint32_t [size_];

	rewind();
	State temp(size_);

	while (sites_>0)
	{
		uncompress(a, b);
		fast_trans(a, c, size_, 32);
		fast_trans(a, d, size_, 32);
		temp.compress(c, d);	
	}
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
	*this=temp;
	cache();	
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
		ret=LZ4_decompress_fast (stream.lz4_ptr, (char *)a, block_size_);
		if (ret > 0)
		{
			stream.lz4_ptr+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		ret=LZ4_decompress_fast (stream.lz4_ptr, (char *)b, block_size_);	
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
	//std::cerr << "started at " << cached_sites_-sites_ << std::endl;
	if ( k < cached_sites_-sites_) rewind();
	while (cached_sites_-sites_ <= k) 
	{
	//	std::cerr << "uncompressed " << cached_sites_-sites_ << std::endl;
		if (sites_==0) {
			fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		int ret=LZ4_decompress_fast (lz4_ptr_, (char *)a, block_size_);

		if (ret > 0)
		{
			lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		ret=LZ4_decompress_fast (lz4_ptr_, (char *)b, block_size_);

		if (ret > 0)
		{
			lz4_ptr_+=ret;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}

		sites_--;
	} 
	//std::cerr << "finished at " << cached_sites_-sites_ << std::endl;
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
State::advance(void)
{
	lz4_ptr_=temp_lz4_ptr_;
	++sites_;
	k_=0;
	cached_=false;
}

void
State::finalize(void)
{
	lz4_ptr_=temp_lz4_ptr_;
	k_=0;
	cached_=false;
}

void 
State::compress_inplace (const uint32_t *a, const uint32_t *b)
{
	int size;
	if (lz4_end_-lz4_ptr_ < 2*block_size_ )
	{
		fprintf(stderr, gettext("mapgd:%s:%d:Increasing buffer size from %f Mb.\n"), __FILE__, __LINE__, float(lz4_end_-lz4_start_)/1000000.);
		increase_buffer_();
	}


	temp_lz4_ptr_=lz4_ptr_;


	size=LZ4_compress_default( (const char*) a, temp_lz4_ptr_, block_size_,  size_t (lz4_end_-temp_lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-temp_lz4_ptr_) );
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	temp_lz4_ptr_+=size;
	size=LZ4_compress_default( (const char*) b, temp_lz4_ptr_, block_size_,  size_t (lz4_end_-temp_lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-temp_lz4_ptr_) );
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	temp_lz4_ptr_+=size;
	if (lz4_end_==temp_lz4_ptr_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}
void 
State::compress (const uint32_t *a, const uint32_t *b)
{
	int size;
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) increase_buffer_();

	size=LZ4_compress_default( (const char*) a, lz4_ptr_, block_size_,  size_t (lz4_end_-lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-lz4_ptr_) );
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	lz4_ptr_+=size;
	size=LZ4_compress_default( (const char*) b, lz4_ptr_, block_size_,  size_t (lz4_end_-lz4_ptr_) > INT_MAX ? INT_MAX : size_t (lz4_end_-lz4_ptr_) );
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

	//std::cerr << "here:" << stream.sites_ << std::endl;
	if (transposed_)
	{
	} else {
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
	}
	delete [] set0;
	delete [] set1;
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
	//std::cerr << "Sites:" << sites_ << std::endl;
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
	delete [] set0;
	delete [] set1;
	cache();
}

std::string 
State::header (void) const
{
	return "@NS:"+std::to_string(size_)+"\tGS:"+std::to_string(sites_)+"\tBS:"+std::to_string(lz4_last_-lz4_start_)+'\n';
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



State sub_sample(const State &state, const size_t &sub_sample, const uint32_t *mask) 
{
	State ret(sub_sample);
	State_stream stream;
	state.set_stream(stream);

	uint32_t *P1=new uint32_t [state.sample_size()];
	uint32_t *P2=new uint32_t [state.sample_size()];
	uint32_t *p1=new uint32_t [sub_sample];
	uint32_t *p2=new uint32_t [sub_sample];

	for (size_t z=0; z<state.genome_size();++z)
	{
		state.uncompress(P1, P2, stream);
		size_t i=0;
		for (size_t x=0; x<state.sample_size(); ++x)
		{
			//std::cerr << x;
			if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			{
				p1[i]=P1[x];
				p2[i++]=P2[x];
			}
			//std::cerr << std::endl;
		}
		ret.compress(p1, p2);
	}
	delete [] P1;
	delete [] P2;
	delete [] p1;
	delete [] p2;
	ret.cache();
	//ret.write(std::cerr);
	return ret;
}
uint8_t
State::get_k(void) const
{
	return k_;
}

void
State::set_k(const uint8_t &k)
{
	k_=k;
};

bool State::empty(void)
{
	return sites_==0;
}
