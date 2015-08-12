/***** interface.h ********************************
 * Description: Accepts arguments, sets flags, autoformats usage, verison, and help
 * Not POSIX compliant...yet.
 * Author: Matthew Ackerman
 * Date: Mon Dec 09 12:00:00 2014.
 * Licence: GNU General Public
 ***************************************************/ 


#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include <list>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <list>
#include <cmath>
#include "interface.h"
#include "../stream_tools.h"
#include "../typedef.h"

int arg_setvectorstr(int, char **, void *);	//<! Takes a reference to a vector<string> and will fill that vector with arguments from arg [] until a string begining with '-' is found.  
int arg_setvectorint(int, char **, void *);	//<! Takes a reference to a vector<int> and calls atoi() on a comma delimited list (no white space) of integers. Ranges of the form N-M are accepted. 
int arg_setint(int, char **, void *);		//<! Takes a reference to an int and calls atoi on the field immediately following the option. Raises errors if the field is not numberic. Returns 2 on success (representing the number of fields proccessed) and -1 on failure. 
int arg_set2int(int, char **, void *);		//<! Takes a reference to a string and sets it to the field immediately following the otpion. Returns 2 on success and -1 on failure. 
int arg_setstr(int, char **, void *);		//<! Takes a reference to ...
int arg_setchar(int, char **, void *);		//<! Same.
int arg_setfloat_t(int, char **, void *);	//<! Same.
int arg_setc_str(int, char **, void *);		//<! Same.

int flag_set(void *);				//<! Sets the reference boolian to true when the option is passed from the command line.
int flag_usage(void *);				//<! Prints an automatically formated usage when called to stdout.
int flag_version(void *);			//<! Prints the version of the environment to stdout. 
int flag_help(void *);				//<! Prints an automatically formated help message to stdout.

int arg_error(int, char **, void *);		//<! Sets the error to be printed if a argument returns -1.
int flag_error(void *);				//<! Sets the error to be printed if a flag_set returns -1.
int command_error(int, char **);		//<! Sets the error to be printed if a command does not return 0.

///breif A class that passes command line flags. It takes no arguments and chews only a single letter.*/
/*! TODO Write a long descrition you jerk!
 */
class flag_t {
	private:
	public:
	flag_t(){
		opt='?';
		lopt="error";
		func=&flag_error;
		emsg="emsg unset";
		umsg="umsg unset";
	};
	flag_t (char opt_, char* lopt_, void *parm_, int (*func_)(void *), char *emsg_, char*umsg_){
		opt=opt_;
		lopt=lopt_;
		func=func_;
		parm=parm_;
		emsg=emsg_;
		umsg=umsg_;
	};
	bool set;			//!< flag toggles whether option has been set.
	char opt;			//!< the option name.
	char *lopt;			//!< the long option name.
	void *parm;			//!< pointer to the parameter to be set.
	int (*func)(void *);		//!< the function to set the parameters.
	char *emsg;			//!< A short error message to display when the proper parameters aren't passed to this option.
	char *umsg;			//!< A short discription of this option to be displayed in the usage message.
};

/// A command line argument. See the interface tutorial for a demonstration of usage. 
/*! TODO Write a long descrition you jerk!
 *
 */
class arg_t {
	private:
	public:
	arg_t(){
		opt='?';
		lopt="error";
		func=&arg_error;
		emsg="emsg unset";
		umsg="umsg unset";
		set=false;
	};
	arg_t(char opt_, char* lopt_, void * parm_, int (*func_)(int, char **, void *), char *emsg_, char*umsg_){
		opt=opt_;
		lopt=lopt_;
		parm=parm_;
		func=func_;
		emsg=emsg_;
		umsg=umsg_;
		set=false;
	};
	bool set;				//!< flag toggles whether option has been set.
	char opt;				//!< the option name.
	char *lopt;				//!< the long option name.
	void *parm;				//!< pointer to the parameter to be set.
	int (*func)(int, char **, void *);	//!< the function to set the parameters.
	char *emsg;				//!< A short error message to display when the proper parameters aren't passed to this option.
	char *umsg;				//!< A short description of this option to be displayed in the usage message.
};

///A sub-command of mapgd.
/*! TODO Write a long descrition you jerk!
 *
 */
class com_t {
	private:
	public:
	com_t(){
		opt='?';
		lopt="error";
		emsg="emsg unset";
		func=&command_error;
		umsg="umsg unset";
		set=false;
	};
	com_t(char opt_, char* lopt_, int (*func_)(int, char **), char *emsg_, char*umsg_){
		opt=opt_;
		lopt=lopt_;
		func=func_;
		emsg=emsg_;
		umsg=umsg_;
		set=false;
	};
	bool set;				/*!< flag toggles whether option has been set*/
	char opt;				/*!< the option name*/
	char *lopt;				/*!< the long option name*/
	int (*func)(int, char **);		/*!< the function to set the parameters*/
	char *emsg;				/*!< A short error message to display when the proper parameters aren't passed to this option*/
	char *umsg;				/*!< A short description of this command to be displated in the usage and help message.*/
};

/// A class for handelling options and automatically formating --help, -h, -u and -v.
/*! TODO Write a long descrition you jerk!
 *
 */
class env_t{
	private:
	public:
	void setname(const char *c)
	{
		name=c;
	};
	void setver(const char *c)
	{
		version=c;
	};
	void setauthor(const char *c){
		author=c;
	};
	void setdescription(const char *c){
		description=c;
	};

	const char *name;				/*!< the name of the command*/
	const char *version;				/*!< the version*/
	const char *author;				/*!< author*/
	const char *description;			/*!< a breif description of the command being executed*/

	std::list <flag_t> flags;			/*!< A list of flags that can be set*/
	std::list <arg_t> args;				/*!< A list of options that can be passed from the command line.*/
	std::list <com_t> commands;			/*!< A list of sub-commands that can be called from the command line.*/

	std::list <arg_t *> required_args;		/*!< A list of flags or options, all of which must be set*/
	std::list <com_t *> required_coms;		/*!< A list of sub-commands one of which must be run*/

	env_t(){
		name="Unnammed program";
		version="0.0";
		author="Unkown";
		description="Unknown purpose";
	};
	
	/*!	\breif adds an optional argument to the list of arguments accepted by the program.*/
	void optional_arg (char opt_, char* lopt_, void * parm_, int (*func_)(int, char **, void *), char *emsg_, char*umsg_)
	{
		args.push_back(arg_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
	};

	/*!	\breif adds an required argument to the list of arguments accepted by the program. -1 is returned by parsargs when required arguments are not provided*/
	void required_arg (char opt_, char* lopt_, void * parm_, int (*func_)(int, char **, void *), char *emsg_, char*umsg_)
	{
		args.push_back(arg_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
		required_args.push_back(&args.back());
	};

	/*!	\breif adds a function to be executed that can accept argc and argv [] arguments. 
	 *
	 *	If a command is set, then an error is raised if no command is provided. Commands execute in order in which they are passed by the user. 
	 * 	Commands only return an exit status, so no further processing of the argc and argv variables can be done by the environment in which a command is called.
	 */
	void command (char opt_, char* lopt_, int (*func_)(int, char **), char *emsg_, char*umsg_)
	{
		commands.push_back(com_t(opt_, lopt_, func_, emsg_, umsg_) );
	};

	/*!	\breif adds a flag to the list of flags that can be accepted by the environment. 
	 *
	 *	Flags can be passed either seperately, or in a concatenated list if thier short form is used.
	 */
	void flag (char opt_, char* lopt_, void * parm_, int (*func_)(void *), char *emsg_, char*umsg_)
	{
		flags.push_back(flag_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
	};

	/*!	\breif Checks to see if all required arguments are set. 
	 *
	 *	Returens 1 if all required arguments are set, 0 if they are not.
	 */
	bool required_set (void)
	{
		std::list <arg_t *>::iterator rarg=required_args.begin();
		std::list <arg_t *>::iterator end=required_args.end();
		while(rarg!=end){
			if(!(*rarg)->set) {std::cerr << (*rarg)->emsg << std::endl; return false;}
			++rarg;
		};
		return true;
	};
	void close(void);
};

/*!	\breif parses the command line arguments. 
 *
 *	Returens 0 if all went well, -1 the rest of the time.
 */
int parsargs(int argc, char *argv[], env_t env);

/*!	\breif prints the automatically generated help file and exist. 
 */
void printHelp(env_t env);

/*!	\breif prints the automatically generated version line. 
 */
void printVersion(env_t env);

/*!	\breif prints the automatically generated usage.
 *	
 *	Currently usage is not correctly formated. Need to go back and make it more helpful.  
 */
void printUsage(env_t env);
#endif 
