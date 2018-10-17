
#include "mapgd_help.h"

//! The help command.
/* mapgd's help command. This may be switched over to simply call up man pages,
 * but for now it prints out formated lines from a flat file.
 */
int mapgd_help(int argc, char *argv[])
{
	if (argc!=2) {
		std::cerr << "No manual entry for ''\nPlease type mapgd help [CONCEPT]\n";
	} else {
#ifdef USE_MAN
		char call[255];
		snprintf(call, 200, "man -M %s mapgd-%s", MANPATH, argv[1]);
		return system(call);
#else
		std::cerr << "mapgd_help is currently unavailable without the command 'man'.\n";
#endif
	}
	
	return 0;					//Since everything worked, return 0!.
}
