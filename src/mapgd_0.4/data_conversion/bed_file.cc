#include "bed_file.h"

#include <bitset>

Bed_data::Bed_data ()
{
	buffer_=NULL;
	size_=0;
}

Bed_data::Bed_data(const size_t &N)
{
	size_=(N >> 2)+((N & 0xC0)!=0);
	std::cerr << ":"<< size_ << std::endl;
	buffer_=new uint8_t[size_];
}

void 
Bed_file::read (Bed_data &bed)
{
	if (open_) 
	{
		if ( (in_->read( (char *)(bed.buffer_), bed.size_) ).eof() ) 
		{
			fprintf(stderr, gettext("mapgd:%s:%d: attempt to read from file failed: eof.\n"), __FILE__, __LINE__ );
			close_table();
		}
		
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: attempt to read from unopened file.\n"), __FILE__, __LINE__ );
		exit(0);
	}
}

#ifndef NOLZ4
id1_t
Bed_data::get(State &state) const
{
	size_t N=state.sample_size();

	uint32_t *s1=new uint32_t [N];
	uint32_t *s2=new uint32_t [N];

	uint32_t *s3=new uint32_t [N];
	uint32_t *s4=new uint32_t [N];


	int k=state.get_k();

	if (k!=32 && k!=0) 
	{
		state.uncompress_inplace(s1, s2, s3, s4);
	} else {
		state.advance();
		memset(s1, 0, sizeof(int32_t)*N);
		memset(s2, 0, sizeof(int32_t)*N);
		memset(s3, 0xff, sizeof(int32_t)*N);
		memset(s4, 0xff, sizeof(int32_t)*N);
		k=0;
	}

	size_t x=0;

	uint8_t *git=buffer_;
	uint8_t *end=buffer_+size_;

	for (; git<end; git++)
	{
		for (int y=0; y<8; y+=2)
		{
			if (x<N){
				switch ( ( *git & (0x03 << y) ) >> y )
				{
					case(0):
					break;
					case(1):
						s3[x] ^= 1 << k;
						s4[x] ^= 1 << k;
					break;
					case(2):
						s1[x] |= 1 << k;
					break;
					case(3):
						s1[x] |= 1 << k;
						s2[x] |= 1 << k;
					break;
					default:
						fprintf(stderr, gettext("mapgd:%s:%d: impossible error.\n"), __FILE__, __LINE__ );
						exit(0);
					break;
				}
			}
			x++;
		}
	}

	state.compress_inplace(s1, s2, s3, s4);
	k++;
	state.set_k(k);

	delete [] s1;
	delete [] s2;
	delete [] s3;
	delete [] s4;

	return 0;
}
#endif

void 
Bed_data::get(Data *data, ...) const
{
	if (0)
	{
	//	get_text("GOOD WORK %s", "BOO!");
	} else {
	//	get_text("BAD WORK %s", "BOO!");
	}
}

void
Bed_data::put(const Data *data, ...)
{
	va_list args;
	va_start(args, data);
}

void
Bed_file::open(const std::ios_base::openmode &mode)
{
//	else if (mode & WRITE) open(std::cout, mode);
//	else {
		fprintf(stderr, gettext("mapgd:%s:%d: cannot open file in std::ios_base::openmode &mode=%d.\n"), __FILE__, __LINE__, mode );
//	}
}


void
Bed_file::open(const char *file_name, const std::ios_base::openmode &mode)
{
	if (mode & READ)
	{
		file_.open(file_name, READ);
		if ( file_) 
		{
			uint32_t bits;
			file_.read((char *)(&bits), 3);
			std::cerr << std::hex <<( bits & 0x00ffffff) << std::endl;

//			if ( (bits & 0x00ffffff)==0x006c1b01) {
			if ( (bits & 0x00ffffff)==0x00011b6c) {
				open_=true;
				table_open_=true;
				in_=&file_;
				return;
			} else {
				fprintf(stderr, gettext("mapgd:%s:%d: Magick bits not set.\n"), __FILE__, __LINE__, file_name );
			}
		}
		
		fprintf(stderr, gettext("mapgd:%s:%d: Could not open file \"%s\" for reading.\n"), __FILE__, __LINE__, file_name );
		exit(EXIT_FAILURE);
	}
	else if (mode & WRITE)
	{
		file_.open(file_name, WRITE);
		if ( file_) open_=true;
	}
	if (open_)
	{
		if (mode & READ)
		{
		}
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Unexpected error opening file \"%s\".\n"), __FILE__, __LINE__, file_name );
	}
		
}
void 
Bed_file::close(void)
{
	if (open_)
	{
		file_.close();
		open_=false;
	}
		
}
