/***** eprintf.h ************************************************************
 * Description: Header file for eprintf, which provides error-
 *              handling capabilities.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Dec 17 11:16:37 2004. 
 * 
 * This file is part of kr.
 *
 * kr is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * kr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with kr; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/
#include <stdio.h>
extern FILE *efopen(char *fname, char *mode); 
extern void eprintf(char *, ...);
extern char *estrdup (char *);
extern void *emalloc(size_t);
extern void *erealloc(void *, size_t);
extern char *progname(void);
extern void setprogname2(char *);
