/***** interface.h ********************************
 * Description: Accepts arguments, sets flags, autoformats usage, verison, and help
 * Not POSIX compliant... yet.
 * Author: Matthew Ackerman
 * Date: Mon Dec 09 12:00:00 2014.
 * Licence: GNU General Public
 ***************************************************/ 


#ifndef INTERFACE_H_
#define INTERFACE_H_	1_
#include <list>
//TODO: this is just here for debuging.
#include <iostream>

int arg_setvectorstr(int, char **, void *);
int arg_setvectorint(int, char **, void *);
int arg_setint(int, char **, void *);
int arg_set2int(int, char **, void *);
int arg_setstr(int, char **, void *);
int arg_setchar(int, char **, void *);
int arg_setfloat_t(int, char **, void *);
int arg_setc_str(int, char **, void *);

int flag_set(void *);
int flag_usage(void *);
int flag_version(void *);
int flag_help(void *);

int arg_error(int, char **, void *);
int flag_error(void *);
int command_error(int, char **);

/*Flag class */
class flag_t {
	//A class for command line flags (options which take no argument) 
	private:
	public:
	flag_t(){
		opt='?';
		lopt="error";
		func=&flag_error;
		emsg="emsg unset";
		umsg="umsg unset";
	};
	flag_t(char opt_, char* lopt_, void *parm_, int (*func_)(void *), char *emsg_, char*umsg_){
		opt=opt_;
		lopt=lopt_;
		func=func_;
		parm=parm_;
		emsg=emsg_;
		umsg=umsg_;
	};
	bool set;			/*flag toggles whether option has been set*/
	char opt;			/*the option name*/
	char *lopt;			/*the long option name*/
	void *parm;			/*pointer to the parameter to be set*/
	int (*func)(void *);		/*the function to set the parameters*/
	char *emsg;			/*A short error message to display when the proper parameters aren't passed to this option*/
	char *umsg;			/*A short discription of this option to be displayed in the usage message*/	
};

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
	bool set;				/*flag toggles whether option has been set*/
	char opt;				/*the option name*/
	char *lopt;				/*the long option name*/
	void *parm;				/*pointer to the parameter to be set*/
	int (*func)(int, char **, void *);	/*the function to set the parameters*/
	char *emsg;				/*A short error message to display when the proper parameters aren't passed to this option*/
	char *umsg;				/*A short discription of this option to be displayed in the usage message*/	
};

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
	bool set;				/*flag toggles whether option has been set*/
	char opt;				/*the option name*/
	char *lopt;				/*the long option name*/
	int (*func)(int, char **);		/*the function to set the parameters*/
	char *emsg;				/*A short error message to display when the proper parameters aren't passed to this option*/
	char *umsg;				/*A short error message to display when the proper parameters aren't passed to this option*/
};

/*information that is displayed in usage and version*/
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
	const char *name;				/*the name of the command*/
	const char *version;				/*the version*/
	const char *author;				/*author*/
	const char *description;			/*a breif description*/

	std::list <flag_t> flags;
	std::list <arg_t> args;
	std::list <com_t> commands;

	std::list <arg_t *> required_args;
	std::list <com_t *> required_coms;

	env_t(){
		name="Unnammed program";
		version="0.0";
		author="Unkown";
		description="Unknown purpose";
	};

	void optional_arg (char opt_, char* lopt_, void * parm_, int (*func_)(int, char **, void *), char *emsg_, char*umsg_)
	{
		args.push_back(arg_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
	};

	void required_arg (char opt_, char* lopt_, void * parm_, int (*func_)(int, char **, void *), char *emsg_, char*umsg_)
	{
		args.push_back(arg_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
		required_args.push_back(&args.back());
	};

	void command (char opt_, char* lopt_, int (*func_)(int, char **), char *emsg_, char*umsg_)
	{
		commands.push_back(com_t(opt_, lopt_, func_, emsg_, umsg_) );
	};

	void flag (char opt_, char* lopt_, void * parm_, int (*func_)(void *), char *emsg_, char*umsg_)
	{
		flags.push_back(flag_t(opt_, lopt_, parm_, func_, emsg_, umsg_) );
	};

	bool required_set (void)
	{
		std::list <arg_t *>::iterator rarg=required_args.begin();
		std::list <arg_t *>::iterator end=required_args.end();
		while(rarg!=end){
			std::cout << "HI!\n";
			if(!(*rarg)->set) {std::cerr << (*rarg)->emsg << std::endl; return false;}
			std::cout << (*rarg)->opt << " is set.\n";
			++rarg;
		};
		return true;
	};
	void close(void);
};

int parsargs(int argc, char *argv[], env_t env);
void printHelp(env_t env);
void printVersion(env_t env);
void printUsage(env_t env);
#endif 
