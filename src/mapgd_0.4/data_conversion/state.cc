#include "state.h"

const std::string State::file_name=".bgt";
const std::string State::table_name="STATE";
const bool State::binary=true;

const Registration State::registered=Registration(State::table_name, State::create);

State::State ()
{
	lz4_start_ = new char [LZ4_BUFFER_SIZE];
	//std::cerr << "lz4_start_ = " << std::hex << size_t(lz4_start_) << std::endl;
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
	if (column_names.size()==3) 
	{
		uint32_t sample_size=atoi(split(column_names[0], ':')[1].c_str() );
		uint32_t sites_size=atoi(split(column_names[1], ':')[1].c_str() );
		uint32_t buffer_size=atoi(split(column_names[2], ':')[1].c_str() );

		size_=sample_size;
		sites_=sites_size;
		cached_sites_=sites_size;

		masked_=MASKED;
		cached_=true;

		lz4_buffer_size_=LZ4_BUFFER_SIZE;

		lz4_start_ = new char [lz4_buffer_size_];
		//std::cerr << "lz4_start_ = " << std::hex << size_t(lz4_start_) << std::endl;
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
	cached_=true;
	masked_=MASKED;

	lz4_buffer_size_=LZ4_BUFFER_SIZE;
	lz4_start_ = new char [lz4_buffer_size_];
	//std::cerr << "lz4_start_ = " << std::hex << size_t(lz4_start_) << std::endl;
	lz4_ptr_=lz4_start_;
	lz4_last_ = lz4_start_ + buffer_size;
	lz4_end_ = lz4_start_ + lz4_buffer_size_;
	temp_lz4_ptr_=lz4_ptr_;

	size_=sample_size;
	sites_=sites_size;
	cached_sites_=sites_size;
	block_size_=sizeof(uint32_t)*size_;

	k_=0;
}

State::State(const uint32_t &sample_size)
{
	transposed_=false;
	masked_=MASKED;
	cached_=false;

	size_=sample_size;
	sites_=0;
	cached_sites_=0;
	lz4_buffer_size_=LZ4_BUFFER_SIZE;

	lz4_start_ = new char [lz4_buffer_size_];
	//std::cerr << "lz4_start_ = " << std::hex << size_t(lz4_start_) << std::endl;
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
		//std::cerr << "Death to " << size_t(lz4_start_) << std::endl;
		clear();
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
		uncompress_(temp_lz4_ptr_, a, b);
		if(masked_) uncompress_(temp_lz4_ptr_);
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
		uncompress_(lz4_ptr_, a, b);
		if (masked_) uncompress_(lz4_ptr_);
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::uncompress (uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d)
{
	if (!cached_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from uncached stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	if (sites_-- ) 
	{
		uncompress_(lz4_ptr_, a, b);
		if (masked_)
		{
			uncompress_(lz4_ptr_, c, d);
		}
		else {
			memset( (char *)c, 0xff, sizeof(uint32_t)*size_); 
			memset( (char *)d, 0xff, sizeof(uint32_t)*size_); 
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
State::uncompress_(char *&ptr, uint32_t *&a, uint32_t *&b) const
{
	int ret=LZ4_decompress_fast (ptr, (char *)a, block_size_);
	if (ret > 0)
	{
		ptr+=ret;
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	ret=LZ4_decompress_fast (ptr, (char *)b, block_size_);	
	if (ret > 0)
	{
		ptr+=ret;
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Malformed LZ4 block. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void 
State::uncompress_(char * &ptr ) const
{
	uint32_t *temp=new uint32_t[size_];
	uncompress_(ptr, temp, temp);
	delete [] temp;
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
		uncompress_(stream.lz4_ptr, a, b);
		if(masked_) uncompress_(stream.lz4_ptr);
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);

	}
}

void 
State::uncompress (uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d, State_stream &stream) const
{
	if (!cached_)
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from uncached stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	if (stream.sites_-- ) 
	{
		uncompress_(stream.lz4_ptr, a, b);
		if(masked_) uncompress_(stream.lz4_ptr, c, d);
		else {
			memset((char *)c, 0xff, sizeof(uint32_t)*size_ );
			memset((char *)d, 0xff, sizeof(uint32_t)*size_ );
		}
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
		exit(LZ4);

	}
}

void 
State::uncompress (uint32_t *a, uint32_t *b, const uint32_t &k)
{
	if ( k < cached_sites_-sites_) rewind();
	while (cached_sites_-sites_ <= k) 
	{
		if (sites_==0) {
			fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		uncompress_(lz4_ptr_, a, b);
		if(masked_) uncompress_(lz4_ptr_);
		sites_--;
	} 
}

void 
State::uncompress (uint32_t *a, uint32_t *b, uint32_t *c, uint32_t *d, const uint32_t &k)
{
	if ( k < cached_sites_-sites_) rewind();
	while (cached_sites_-sites_ <= k) 
	{
		if (sites_==0) {
			fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from empty stream. Exiting.\n"), __FILE__, __LINE__);
			exit(LZ4);
		}
		uncompress_(lz4_ptr_, a, b);
		if (masked_) uncompress_(lz4_ptr_, c, d);
		else 
		{
			memset((char *)c, 0xff, sizeof(uint32_t)*size_ );
			memset((char *)d, 0xff, sizeof(uint32_t)*size_ );
		}
		sites_--;
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
	std::cerr << "allocating " << lz4_buffer_size_+LZ4_BUFFER_SIZE << " at "  << size_t(new_lz4) << std::endl;
	if (!new_lz4) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Failed to allocate memory. Mother fucker. Exiting.\n"), __FILE__, __LINE__);
	}
	memcpy(new_lz4, lz4_start_, lz4_buffer_size_);
	lz4_buffer_size_+=LZ4_BUFFER_SIZE;
	std::cerr << size_t (lz4_ptr_-lz4_start_) << std::endl;
	std::cerr << size_t (lz4_last_-lz4_start_) << std::endl;
	std::cerr << size_t (lz4_buffer_size_) << std::endl;

	lz4_ptr_=new_lz4+(lz4_ptr_-lz4_start_);
	lz4_last_=new_lz4+(lz4_last_-lz4_start_);
	lz4_end_=new_lz4+lz4_buffer_size_;
	delete [] lz4_start_;
	lz4_start_=new_lz4;
	std::cerr << size_t(lz4_end_-lz4_start_) << '\t' << size_t(lz4_ptr_-lz4_start_) << std::endl;
	std::cerr << size_t(lz4_end_-lz4_ptr_) << '\t' << size_t() << std::endl;
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
State::compress_(char *&ptr, char *&end, const uint32_t *a, const uint32_t *b) const
{
	int size=LZ4_compress_default( (const char*) a, ptr, block_size_,  size_t (end-ptr) > INT_MAX ? INT_MAX : size_t (end-ptr) );
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	ptr+=size;
	size=LZ4_compress_default( (const char*) b, ptr, block_size_,  size_t (end-ptr) > INT_MAX ? INT_MAX : size_t (end-ptr) );
	if (size==0) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
	ptr+=size;
	if (end==ptr)
	{
		fprintf(stderr, gettext("mapgd:%s:%d:FLAGRANT SYSTEM ERROR. Computer over. lz4buffer = Full. (Write me an e-mail!)\n"), __FILE__, __LINE__);
		exit(LZ4);
	}
}

void
State::compress_(char *&ptr, char *&end) const
{
	char * a=new char [block_size_];
	memset(a, 0xff, block_size_);
	compress_(ptr, end, (uint32_t *)a, (uint32_t *)a);
	delete [] a;
}
void 
State::compress_inplace (const uint32_t *a, const uint32_t *b, const uint32_t *c, const uint32_t *d)
{
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) {
		std::cerr << __LINE__ << " compressing in place is increasing buffer...\n";
		increase_buffer_();
	}

	char *temp_lz4_ptr=lz4_ptr_;
	compress_(temp_lz4_ptr, lz4_end_, a, b);
	if(masked_) compress_(temp_lz4_ptr, lz4_end_, c, d);
}

void 
State::compress_inplace (const uint32_t *a, const uint32_t *b)
{
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) 
	{
		std::cerr << __LINE__ << " compress in place is increasing buffer...\n";
		increase_buffer_();
	}

	char *temp_lz4_ptr=lz4_ptr_;
	compress_(temp_lz4_ptr, lz4_end_, a, b);
	if(masked_) compress_(temp_lz4_ptr, lz4_end_);
}

void 
State::compress (const uint32_t *a, const uint32_t *b, const uint32_t *c, const uint32_t *d)
{
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) 
	{
		std::cerr << __LINE__ << " compress is increasing buffer...\n";
		increase_buffer_();
	}

	compress_(lz4_ptr_, lz4_end_, a, b);
	if(masked_) compress_(lz4_ptr_, lz4_end_, c, d);

	++sites_;
	cached_=false;
}

void 
State::compress (const uint32_t *a, const uint32_t *b)
{
	if (lz4_end_-lz4_ptr_ < 2*block_size_ ) 
	{
		
		std::cerr << __LINE__ << " compress is increasing buffer because " << size_t (lz4_ptr_) << " is "  << size_t(lz4_end_-lz4_ptr_) << '\t' << 2*block_size_ << std::endl;
		increase_buffer_();
	}

	compress_(lz4_ptr_, lz4_end_, a, b);
	if(masked_) compress_(lz4_ptr_, lz4_end_);

	++sites_;
	cached_=false;
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
				int c=0;
				for (size_t y=0; y<size_; ++y)
				{
					//std::cerr  << std::string(c, '\t') << ( (set0[y] & mask[b]) >> b) <<  ( ( set1[y] & mask[b] ) >>b );
					out  << std::string(c, '\t') << ( (set0[y] & mask[b]) >> b) <<  ( ( set1[y] & mask[b] ) >>b );
					/*if (set0[y] & mask[b]) {
						if ( (set1[y] & mask[b] ) ) 
						{
							out << std::string(c, '\t') <<"1	1";
							} else {
							out << std::string(c, '\t') <<"1	0";
						}
					} else {
						if ( (set1[y] & mask[b]) )
						{
							out << std::string(c, '\t') <<"0	1";
						} else {
							out << std::string(c, '\t') <<"0	0";
							}
					}*/
					c=1;
				}
				if (stream.sites_!=0 || b!=31) out << std::endl;
			}
		}
	}
	delete [] set0;
	delete [] set1;
}

/*
void 
State::write (std::ostream& out) const 
{
	uint32_t *set0=new uint32_t[size_];
	uint32_t *set1=new uint32_t[size_];

	uint32_t *set2=new uint32_t[size_];
	uint32_t *set3=new uint32_t[size_];

//	memset((char *)set2, 0xff, sizeof(uint32_t)*size_ );
//	memset((char *)set3, 0xff, sizeof(uint32_t)*size_ );

	State_stream stream;
	set_stream(stream);

	char o[]="o1i0";
	//std::cerr << "here:" << stream.sites_ << std::endl;
	if (transposed_)
	{
	} else {
		while (stream.sites_>0)
		{
			if(masked_) uncompress(set0, set1, set2, set3, stream);
			else uncompress(set0, set1, stream);
			for (size_t b=0; b<32; ++b) 
			{
				int c=0;
				for (size_t y=0; y<size_; ++y)
				{
					short d=(set0[y] & mask[b] >> b) << 1;
					short e=(set1[y] & mask[b] >> b) << 1;
					short f=(set2[y] & mask[b] >> b);
					short g=(set3[y] & mask[b] >> b);
					//std::cerr << d << ", " << e << ", " << f << ", " << g << std::endl;
					out << std::string(c, '\t') << o[d+f] << o[e+g];
					c=1;
				}
				if (stream.sites_!=0 || b!=31) out << std::endl;
			}
		}
	}
	delete [] set0;
	delete [] set1;

	delete [] set2;
	delete [] set3;
}*/
	
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
	uint32_t *set2=new uint32_t[size_];
	uint32_t *set3=new uint32_t[size_];

	char w[2], z[2]; 
	uint32_t sites=sites_;
	//std::cerr << "Sites:" << sites_ << std::endl;
	clear();
	for (size_t x=0; x<sites; ++x)
	{
		memset(set0, 0, sizeof(uint32_t)*size_);
		memset(set1, 0, sizeof(uint32_t)*size_);
		memset(set2, 0xff, sizeof(uint32_t)*size_);
		memset(set3, 0xff, sizeof(uint32_t)*size_);


		for (size_t b=0; b<32; ++b) 
		{
			for (size_t y=0; y<size_; ++y)
			{
				in >> w;
				if (w[0]=='i' || w[0]=='1' ) set0[y]+=mask[b];
				if (w[1]=='i' || w[1]=='1' ) set1[y]+=mask[b];
				if (w[0]=='i' || w[0]=='l' ) set0[y]-=mask[b];
				if (w[1]=='i' || w[1]=='l' ) set1[y]-=mask[b];
			}
		}
		compress(set0, set1, set2, set3);
	}
	delete [] set0;
	delete [] set1;

	delete [] set2;
	delete [] set3;

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



State 
sub_sample(const State &state, const size_t &sub_sample, const uint32_t *mask) 
{
	State ret(sub_sample);
	State_stream stream;
	state.set_stream(stream);

	uint32_t *P1=new uint32_t [state.sample_size()];
	uint32_t *P2=new uint32_t [state.sample_size()];
	uint32_t *p1=new uint32_t [sub_sample];
	uint32_t *p2=new uint32_t [sub_sample];

	uint32_t *M1=new uint32_t [state.sample_size()];
	uint32_t *M2=new uint32_t [state.sample_size()];
	uint32_t *m1=new uint32_t [sub_sample];
	uint32_t *m2=new uint32_t [sub_sample];

	for (size_t z=0; z<state.genome_size();++z)
	{
		state.uncompress(P1, P2, M1, M2, stream);
		size_t i=0;
		for (size_t x=0; x<state.sample_size(); ++x)
		{
			if ( mask[x >> 5] & (1 << (x & 0x1F) ) ) 
			{
				p1[i]=P1[x];
				p2[i]=P2[x];
				m1[i]=M1[x];
				m2[i++]=M2[x];
			}
		}
		ret.compress(p1, p2, m1, m2);
	}

	delete [] P1;
	delete [] P2;
	delete [] p1;
	delete [] p2;

	delete [] M1;
	delete [] M2;
	delete [] m1;
	delete [] m2;

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
}

bool State::empty(void)
{
	return sites_==0;
}
