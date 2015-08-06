#include "basic_types.h"

char* CHAR7::to_text(uint8_t *source)
{
	char text[8];
	memcpy(text, source, 7);
	memset(text+7, 0);
	return text;
}

void CHAR7::from_text(char *source, uint8_t *destination)
{
	memcpy(destination, source, 7);
	memset(destination+7, 0);
	return text;
}

char* CHAR13::to_text(uint8_t *source)
{
	char text[14];
	memcpy(text, source, 13);
	memset(text+13, 0);
	return text;
}

void CHAR13::from_text(char *source, uint8_t *destination)
{
	memcpy(destination, source, 13);
	memset(destination+13, 0);
	return text;
}

char* UINT8::to_text(uint8_t *source)
{
	char text[8];
	size_t wrote=snprintf(text, 8, "%d", *source);
	if (wrote<=7) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": warning: number too large for 7 character representation.";
	}
	memset(text+7, 0);
	return text;
}

void UINT8::from_text(char *source, uint8_t *destination)
{
	char text[8];
	memcpy(text, source, 7);
	memset(text+7, 0);
	int value=atoi(text);
	if (value<0 || value>255) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": warning: out of range.";
	}
	uint8_t ret=value;
	return ret;
}
