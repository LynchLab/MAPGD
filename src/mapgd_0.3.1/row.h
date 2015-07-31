///A class to read store abstracted data in an ... 
class row{
private:
	void *columns_;
	key_t *column_keys_;
	size_t size_binary_;
	size_t size_text_;
	char delim_;
	read_row_text();	
	read_row_binary();	
	write_row_text();	
	write_row_binary();	
public:
	row 
	void* get(const key_t &);
	void put(const key_t &, void*);
};

read_row(istream &in, row &);
write_row(osteam &out, const row &);
