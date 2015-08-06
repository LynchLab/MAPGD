#include "map-file.h"

bool is_open(void) const
{
	return open_;
}

void map_seekg(id1_t pos)
{
	file.seekg((streampos)(pos*size_), std::ios_base::beg);
}

void map_seekp(id1_t pos)
{
	file.seekp((streampos)(pos*size_), std::ios_base::beg);
}

void map_seekg(id1_signed_t off, std::ios_base::seekdir dir)
{
	file.seekg((streampos)(off*size_), dir);
}

void map_seekp(id1_signed_t off, std::ios_base::seekdir dir)
{
	file.seekp((streampos)(off*size_), dir);
}

id1_t tellp(void)
{
	return pos_p_;
}

id1_t tellg(void)
{
	return pos_g_;
}

map_file* open(const char *, std::ios_base::openmode)
{
}
map_file* open(std::ios_base::openmode)
{
}
void close(void)
{
}

key get_key(char *key_name)
{
}
std::list <key> get_keys(void)
{
}

void format(std::list <key>);                           //!< formats the file. This can only be done if row is not initalized.

        /*functions dealing with ?*/
        size_t size(void) const;                                //!< Retuern the number of rows in the file. Returns 0 if unknown.
                                                                // If pro file is opened from a std::in, then the number of rows is simply 
                                                                // the number of rows that have been read in until d. However, if 

        size_t column_size(void) const;                         //!< Returns the number of columns in the file. Returns 0 if unknown.

        file_index get_index(void) const;               //!< Returns the file_index.
        const id1_t get_line_number(void) const;        //!< 


