#include "basic_types.h"

void CHAR7::to_text(uint8_t *source, char *destination)
{
	memcpy(destination, source, 7);
	memset(destination+7, 0, 1);
}

void CHAR7::from_text(char *source, uint8_t *destination)
{
	memcpy(destination, source, 7);
	memset(destination+7, 0, 1);
}

void CHAR13::to_text(uint8_t *source, char *destination)
{
	memcpy(destination, source, 13);
	memset(destination+13, 0, 1);
}

void CHAR13::from_text(char *source, uint8_t *destination)
{
	memcpy(destination, source, 13);
	memset(destination+13, 0, 1);
}

void UINT8::to_text(uint8_t *source, char *destination)
{
	size_t wrote=snprintf(destination, 7, "%d", *source);
	if (wrote<=7) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": warning: number too large for 7 character representation.";
	}
	
}

void UINT8::from_text(char *source, uint8_t *destination)
{
	char buffer[8]={0};
	memcpy(buffer, source, 7);
	int value=atoi(buffer);
	if (value<0 || value>255) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": warning: out of range.";
	}
	*destination=value;
}
