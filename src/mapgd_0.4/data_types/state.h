#ifndef _STATE_H_
#define _STATE_H_

#include <string.h>
#include <iostream>
#include <sstream>
#include <climits>
#include <libintl.h>

#include "lz4.h"		//Fast compression/decompression.

#include "data.h"
#include "typedef.h"
#include "gzstream.h"		//Good compression.
#include "stream_tools.h"

//			964799

#define LZ4_BUFFER_SIZE	2000000000

class State_stream 
{
public:
	char *lz4_ptr;
	size_t size_, sites_, cached_sites_;
};


class State : public virtual Data
{
private:
	size_t size_, sites_, cached_sites_, block_size_, lz4_buffer_size_;
	uint8_t k_;

	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new State(Columns);
	}

	bool cached_, transposed_;
	void increase_buffer_(void);
	
	char *lz4_ptr_, *lz4_start_, *lz4_end_, *lz4_last_, *temp_lz4_ptr_;  
public:
	State();		//!< constructor needed by map_file. String should be column names. 
	State(const std::vector <std::string> &);	

	State(const State &rhs)
	{
		std::cerr << "Copy constructor called from:" << size_t(&rhs) << " to " << size_t (this) << " lz4_buffer_size: " <<  lz4_buffer_size_ << std::endl;
		lz4_buffer_size_=rhs.lz4_buffer_size_;

 		lz4_start_ = new char [lz4_buffer_size_];
		std::cerr << "Allocating " << size_t(lz4_start_) << ", " << lz4_buffer_size_ << std::endl;
		lz4_end_ = lz4_start_ + lz4_buffer_size_;
 		lz4_last_ = lz4_start_+(rhs.lz4_last_-rhs.lz4_start_);
		memcpy(lz4_start_, rhs.lz4_start_, rhs.lz4_end_-rhs.lz4_start_);

		cached_=rhs.cached_;
		transposed_=rhs.transposed_;
		block_size_=rhs.block_size_;
		size_=rhs.size_;
		sites_=rhs.sites_;
		cached_sites_=rhs.cached_sites_;
		
	};

	State(State&& rhs)
	{
		std::cerr << "Move constructor called " << std::endl;
		if (this != &rhs)
		{
			lz4_buffer_size_=rhs.lz4_buffer_size_;
 			lz4_start_ = rhs.lz4_start_;
			lz4_end_ = rhs.lz4_end_;
 			lz4_last_ = rhs.lz4_last_;
			cached_=rhs.cached_;
			block_size_=rhs.block_size_;
			size_=rhs.size_;
			sites_=rhs.sites_;
			transposed_=rhs.transposed_;
			cached_sites_=rhs.cached_sites_;
			rhs.lz4_start_ = nullptr;
		}
	}

	State& operator=(const State& rhs )
	{ 
		std::cerr << "Assignment constructor called " << std::endl;
		if (this != &rhs)
		{
			if (lz4_start_) delete [] lz4_start_;
			std::cerr << "Assign from:" << size_t(&rhs) << " to " << size_t (this) << " lz4_buffer_size: " <<  lz4_buffer_size_ << std::endl;
			lz4_buffer_size_=rhs.lz4_buffer_size_;
	
	 		lz4_start_ = new char [lz4_buffer_size_];
			std::cerr << "Allocating " << size_t(lz4_start_) << ", " << lz4_buffer_size_ << ", " << size_t(rhs.lz4_last_-rhs.lz4_start_) << ", " << size_t(rhs.lz4_end_-rhs.lz4_start_) << std::endl;
			lz4_end_ = lz4_start_ + lz4_buffer_size_;
 			lz4_last_ = lz4_start_+(rhs.lz4_last_-rhs.lz4_start_);
			memcpy(lz4_start_, rhs.lz4_start_, rhs.lz4_end_-rhs.lz4_start_);

			cached_=rhs.cached_;
			transposed_=rhs.transposed_;
			block_size_=rhs.block_size_;
			size_=rhs.size_;
			sites_=rhs.sites_;
			cached_sites_=rhs.cached_sites_;
		}
		return *this;
	}

	State(const uint32_t &);
	State(const uint32_t &, const uint32_t &, const uint32_t &);
	~State();

	void set_k(const uint8_t &);
	uint8_t get_k (void) const;

	void uncompress (uint32_t *a, uint32_t *b);
	void uncompress_inplace (uint32_t *a, uint32_t *b);
	void uncompress (uint32_t *a, uint32_t *b, const uint32_t &k);
	void uncompress (uint32_t *a, uint32_t *b, State_stream &) const;

	void cache (void);
	void rewind (void);
	void advance (void);
	void finalize (void);

	void transpose (void);
	
	void compress (const uint32_t *a, const uint32_t *b);
	void compress_inplace (const uint32_t *a, const uint32_t *b);

	void clear(void);

	static const std::string file_name;				//!< The dafualt extention for files.
	static const std::string table_name;				//!< Destination table in Db.

	const std::string get_file_name() const;
	const std::string get_table_name() const;

	static const bool binary;				//!< Destination table in Db.

	inline const size_t & sample_size(void) const {return size_;};
	inline const size_t & genome_size(void) const {return sites_;};
	inline const size_t cached(void) const {return sites_;};

	void read(std::istream& str);
	void write(std::ostream& str) const;

	void read_binary(std::istream& str);
	void write_binary(std::ostream& str) const;

	std::string header() const;

	void set_stream(State_stream &) const;

	inline size_t get_free(void) const {return (size_t)(lz4_end_-lz4_last_); };

	double compression_ratio(void) const;

	inline void swap(State &rhs)
	{
		std::swap(size_, rhs.size_);
		std::swap(sites_, rhs.sites_);
		std::swap(cached_sites_, rhs.cached_sites_);
		std::swap(block_size_, rhs.block_size_);
		std::swap(cached_, rhs.cached_);
		std::swap(lz4_ptr_, rhs.lz4_ptr_);
		std::swap(lz4_start_, rhs.lz4_start_);
		std::swap(lz4_end_, rhs.lz4_end_); 
		std::swap(lz4_last_, rhs.lz4_last_);
	}

	bool empty(void);
};
	
State sub_sample(const State &, const size_t &, const uint32_t *mask);

/*
inline void 
swap(State &lhs, State &rhs)
{
	std::swap(lhs.size_, rhs.size_);
	std::swap(lhs.sites_, rhs.sites_);
	std::swap(lhs.cached_sites_, rhs.cached_sites_);
	std::swap(lhs.block_size_, rhs.block_size_);
	std::swap(lhs.cached_, rhs.cached_);
	std::swap(lhs.lz4_ptr_, rhs.lz4_ptr_);
	std::swap(lhs.lz4_start_, rhs.lz4_start_);
	std::swap(lhs.lz4_end_, rhs.lz4_end_); 
	std::swap(lhs.lz4_last_, rhs.lz4_last_);
}*/
#endif
