/***** sam2pro.h **********************************
 * Description: Convert sam output to profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 21 22:46:11 2010
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "interface.hpp"
#include "stringUtil.hpp"
#include "eprintf.hpp"
#include "tab.hpp"

int sam2pro(int argc, char *argv[]);
void runAnalysis(FILE *fp, Args *args, int *dic);
void scanCol(char *column, int *dic, int *count, char consensus, char *number);
