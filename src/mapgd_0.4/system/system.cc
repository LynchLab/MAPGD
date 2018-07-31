#include "system.h"

//#ifdef UNIX
size_t system_memory(void)
{
    size_t pages = sysconf(_SC_PHYS_PAGES);
    size_t page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
//#endif 
