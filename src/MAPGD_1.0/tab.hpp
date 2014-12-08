/***** tab.h ******************************************************
 * Description: Interface for tab library that handles tabular
 *              data. In contrast to the original description
 *              (c.f. ref.), this implementation can handle
 *              variable field separators.
 * Reference: Kernighan, B. W. and Pike, R. (1999). The Practice
 *            of programming. Addision Wesley; chapter 4.
 * Author: Bernhard Haubold, bernhard.haubold@fh-weihenstephan.de
 * File created on Wed Dec 22 09:31:39 2004.
 *****************************************************************/

extern char *tabGetLine(FILE *f);                 /* read next input line */
extern char *tabGetPlainLine(FILE *f);            /* get line without splitting */
extern char *tabField(int n);                     /* return field n */
extern int tabNfield(void);                       /* return number of fields */
extern void tabSetFieldSep(char *fieldsep);       /* set the field separator, default: \t */
