/* A command to ensure all table keys are defined in keys.txt  */

#ifndef TEST_KEYS_H_
#define TEST_KEYS_H_

#include <string>
#include <iostream>
#include <sstream>
#include "datatypes.h"
#include "map_file.h"
#include "interface.h"
#include "typedef.h"

/** 
  * \ingroup COMMANDS
  * @{
  */
/// Test all data_types to ensure that values are defined in keys.txt
int test_keys(int argc, char *argv[]);
/** @} */

#endif
