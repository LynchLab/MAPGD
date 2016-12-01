/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line, inspired by a program by Bernard
 * Author: Matthew Ackerman, matthew.s.ackerman@gmail.com
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <cmath>

#include "interface.h"
#include "stream-tools.h"

/* @Breif : A command parser that uses functors. 
* 
*  Tis' a silly thing.
*/


/* @Breif : tells you whether a string is an int */
bool 
isfloat(const char *c)
{
	const char *s=c;
	while (*c !=0 && std::isdigit(*c) ) ++c;
	if (*c !=0 && *c=='.') ++c;
	while (*c !=0 && std::isdigit(*c) ) ++c;
	return (*c==0 && c!=s);
}

/* @Breif : tells you whether a string is an int */
bool 
isint(const char *c)
{
	const char *s=c;
	while (*c !=0 && std::isdigit(*c)) ++c;
	return (*c==0 && c!=s);
}

//! Makes sure that the usage message isn't longer than 79 char.
const char *
format_usage(const char *message, const size_t &padding)
{
	std::string unformated_line(message);
	std::string formated_lines;
	std::string::iterator word;
	while(unformated_line.size()>(79-padding) ){
		word=unformated_line.begin()+(79-padding);
		while (word!=unformated_line.begin() && *word!=' ') word--;
		formated_lines.insert(formated_lines.end(), unformated_line.begin(), word);
		formated_lines.insert(formated_lines.end(), 1, '\n');
		formated_lines.insert(formated_lines.end(), padding, ' ');
		unformated_line.erase(unformated_line.begin(), word+1);
	}
	formated_lines.insert(formated_lines.end(), unformated_line.begin(), unformated_line.end() );
	return formated_lines.c_str();
}

/*@Brief : sets a vector of strings from a string*/
int 
arg_set_vector_str(int argc, char **argv, void *parm)
{
	int n=0;
	std::vector <std::string> *v=(std::vector <std::string> *)(parm);
	while (n<argc){
		if (argv[n][0]=='-') return n;
		v->push_back(argv[n]);
		if (n+1==argc) return n+1;
		++n;
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << ": error parsing " << argv[1] << std::endl;
	exit(1);
}

int 
arg_set_region(int argc, char **argv, void *parm)
{
	if (argc>0){
		if (isregion(argv[0])){
			*( (Region *)(parm) )=ator(argv[0]);
			return 1;
		}
		return ARG_ERROR; 
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << " arg_set_region:error parsing " << argv[1] << std::endl;
	exit(1);
}

/*@Breif : sets a vector of ints from a string. */
int
arg_set_vector_ui(int argc, char **argv, void *parm)
{
	std::vector <unsigned int> *v=(std::vector <unsigned int> *)(parm);
	if (argc>0){
		std::vector<std::string> elems=split(argv[0], ',');
		for (size_t x=0; x<elems.size(); ++x){
			if ( isint(elems[x].c_str() ) ) v->push_back(atoi(elems[x].c_str() )-1) ;
			else {
				std::vector<std::string> intpair=split(elems[x], '-');
				if (intpair.size()==2){
					if (isint(intpair[0].c_str() ) && isint(intpair[1].c_str() ) ){
						for (unsigned int y=atoi(intpair[0].c_str() ); y<=atoi(intpair[1].c_str() ); ++y){
							if (std::find(v->begin(), v->end(), y-1)==v->end() ){
								v->push_back(y-1);
							}
						};
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					};
				} else {
					std::vector<std::string> intpair=split(elems[x], '*');
					if (intpair.size()==2){
						for (unsigned int y=0; y<atoi(intpair[1].c_str() ); ++y){
							v->push_back(atoi(intpair[0].c_str() ) );
						}
						
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					}	
				};
			}
		}
		return 2;
	} 
	std::cerr << "arg_set_vector_uint32: error parsing " << argv[0] << std::endl;
	exit(1);
}

int
arg_set_vector_uli(int argc, char **argv, void *parm)
{
	std::vector <unsigned long int> *v=(std::vector <unsigned long int> *)(parm);
	if (argc>0){
		std::vector<std::string> elems=split(argv[0], ',');
		for (size_t x=0; x<elems.size(); ++x){
			if ( isint(elems[x].c_str() ) ) v->push_back(atol(elems[x].c_str() )-1);
			else {
				std::vector<std::string> intpair=split(elems[x], '-');
				if (intpair.size()==2){
					if (isint(intpair[0].c_str() ) && isint(intpair[1].c_str() ) ){
						for (unsigned long int y=atol(intpair[0].c_str() ); y<=atol(intpair[1].c_str() ); ++y){
							if (std::find(v->begin(), v->end(), y-1)==v->end() ){
								v->push_back(y-1);
							}
						};
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					};
				} else {
					std::vector<std::string> intpair=split(elems[x], '*');
					if (intpair.size()==2){
						for (unsigned int y=0; y<atoi(intpair[1].c_str() ); ++y){
							v->push_back(atoi(intpair[0].c_str() ) );
						}
						
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					}	
				}
			}
		}
		return 1;
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << " arg_set_vector_uint64_t: error parsing " << argv[1] << std::endl;
	exit(1);
}

int
arg_set_vector_ulli(int argc, char **argv, void *parm)
{
	std::vector <unsigned long long int> *v=(std::vector <unsigned long long int> *)(parm);
	if (argc>0){
		std::vector<std::string> elems=split(argv[0], ',');
		for (size_t x=0; x<elems.size(); ++x){
			if ( isint(elems[x].c_str() ) ) v->push_back(atol(elems[x].c_str() )-1);
			else {
				std::vector<std::string> intpair=split(elems[x], '-');
				if (intpair.size()==2){
					if (isint(intpair[0].c_str() ) && isint(intpair[1].c_str() ) ){
						for (unsigned long long int y=atol(intpair[0].c_str() ); y<=atol(intpair[1].c_str() ); ++y){
							if (std::find(v->begin(), v->end(), y-1)==v->end() ){
								v->push_back(y-1);
							}
						};
					} else {
						std::cerr << __FILE__ << ":" << __LINE__ << " cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					};
				} else {
					std::vector<std::string> intpair=split(elems[x], '*');
					if (intpair.size()==2){
						for (unsigned int y=0; y<atoi(intpair[1].c_str() ); ++y){
							v->push_back(atoi(intpair[0].c_str() ) );
						}
						
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					}	
				};
			}
		}
		return 1;
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << " arg_set_vector_uint64_t: error parsing " << argv[1] << std::endl;
	exit(1);
}

// A functor to set integer parameters
int
arg_set_int(int argc, char **argv, void *parm)
{
	if (argc>0){
		if (isint(argv[0])){
			*( (int *)(parm) )=atoi(argv[0]);
			return 1;
		}
		return ARG_ERROR; 
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << " arg_setint:error parsing " << argv[1] << std::endl;
	exit(1);
}

//A functor to set single character parameters
int
arg_set_char(int argc, char **argv, void *parm)
{
	if (argc>0){
		*(char *)(parm)=argv[0][0];
		return 1;
	} exit(1);	
}

//A functor to set strings. Returns the number of fields processed on success (1) and exist on failure.
int
arg_set_str(int argc, char **argv, void *parm)
{
        if (argc>0){
                *(std::string *)(parm)=argv[0];
                return 1;
        } exit(1);
}

//A functor to set c strings. Returns the number of fields processed on success (1) and exist on failure.
int 
arg_setc_str(int argc, char **argv, void *parm)
{
	if (argc>0){
		*(char **)(parm)=argv[0];
		return 1;
	} exit(1);	
}

//A functor to set floats. Returns the number of fields processed on success (1), and exist on failure.
int 
arg_set_float(int argc, char **argv, void *parm)
{
	if (argc>0){
		if (isfloat(argv[0])){
			*(float *)(parm)=atof(argv[0]);
			return 1;
		}
		return ARG_ERROR; 
	} exit(1);	
}

//A functor to set doubles. Returns the number of fields processed on success (1), and exits on failure.
int 
arg_set_double(int argc, char **argv, void *parm)
{
	if (argc>0){
		if (isfloat(argv[0])){
			*(double *)(parm)=atof(argv[0]);
			return 1;
		}
		return ARG_ERROR; 
	} exit(1);	
}

//A functor to set long doubles.
int
arg_set_long_double(int argc, char **argv, void *parm)
{
	if (argc>0){
		if (isfloat(argv[0])){
			*(long double *)(parm)=atof(argv[0]);
			return 1;
		}
		return ARG_ERROR; 
	} exit(1);	
}

//TODO Comment
int 
flag_set(void *parm)
{
	*(bool *)(parm)=!(*(bool *)(parm));
	return 0;
}

//TODO Comment
int 
flag_version(void *parm)
{
	print_version(*(Environment *) parm);
	return 0;
}


//! prints the long help statment to stdout (print_help).
int
flag_help(void *parm)
{
	print_help(*(Environment *) parm);
	return 0;
}

//TODO Comment
int
flag_usage(void *parm)
{
	print_usage(*(Environment *) parm);
	return 0;
}

//TODO Comment
int
arg_error(int argc, char **argv, void *parm)
{
	std::cerr << __FILE__ << ":" << __LINE__ << " error: option functor unset\n";
	return 1;
}

/* The main command line argument parsing routine */

int 
parsargs(int argc, char *argv[], Environment &env)
{
	int optind=1,  arg_ret=0; 
	char *optopt; 

	std::list <Argument>::iterator arg=env.args.begin();
	std::list <Argument>::iterator arg_end=env.args.end();

	std::list <Argument *>::iterator pos=env.positional_args.begin();
	std::list <Argument *>::iterator pos_end=env.positional_args.end();

	std::list <Command>::iterator com=env.commands.begin();
	std::list <Command>::iterator com_end=env.commands.end();

	std::list <Flag>::iterator flag=env.flags.begin();
	std::list <Flag>::iterator flag_end=env.flags.end();

	while (optind<argc){
		if (argv[optind][0]=='-'){
			if (argv[optind][1]!=0){
				if (argv[optind][1]!='-'){
					optopt=argv[optind]+1;
					arg=env.args.begin();
					while(arg!=arg_end){
						if (*optopt==arg->opt){
							if (argv[optind][2]!=0){
								//TODO CHECK THE SUPER DANGERS SHIT HERE.
								char* to=*(argv+optind);
								char* from=*(argv+optind)+2;
								do {*to=*from; ++to, ++from;} while((*from)!=0);
								*to=0;
								arg_ret=arg->func(argc-optind, argv+optind, arg->parm);
							} else {
								optind++;
								arg_ret=arg->func(argc-optind, argv+optind, arg->parm);
							}
							if(arg_ret!=ARG_ERROR) {
								optind+=arg_ret;
								arg->set=true;
							} else {
								std::cerr << env.name << ":" << " option --" << arg->lopt << " -" << arg->opt << " "<< arg->emsg  << std::endl; 
							}
							//TODO insert error evaluation...
							break;
						} ++arg;
					} if(arg==arg_end) {
						while(*optopt!='\0'){
							flag=env.flags.begin();
							while(flag!=flag_end){
								if (*optopt==flag->opt){
									flag->func(flag->parm);
									optopt++;
									break;
								};
								++flag;
							};
							if (*optopt!='\0'){
								if(flag==flag_end) {
									std::cerr << env.name << ":" << " invalid option -- '" << optopt << "'" << std::endl; 
									return -1; 
								};
							}
						}
						optind++;
					} 
				} else {
					optopt=argv[optind]+2;
					arg=env.args.begin();
					while(arg!=arg_end){
						if (strcmp(optopt, arg->lopt)==0){
							optind+=arg->func(argc-optind, argv+optind, arg->parm);
							arg->set=true;
							break;
						} ++arg;
					} if(arg==arg_end){
						flag=env.flags.begin();
						while(flag!=flag_end){
							if (strcmp(optopt, flag->lopt)==0){
								optind+=flag->func(flag->parm);
							};
							++flag;
						};
					} if(flag==flag_end) {
						std::cerr << env.name << ":" << " invalid option -- '" << optopt << "'" << std::endl; 
						return -1; 
					}; 
				}
			} else {
				std::cerr << env.name << ":" << " invalid option -- ' ' " << std::endl; 
				return -1; 
			}; 
		}
		else {
			optopt=argv[optind];
			com=env.commands.begin();
			while(com!=com_end){
				if (strcmp(optopt, com->lopt)==0){
					return com->func(argc-optind, argv+optind);
				}
				++com;
			};
			if (com==com_end) {
				if (pos!=pos_end){
					optind+=(*pos)->func(argc-optind, argv+optind, (*pos)->parm);
					optind++;
					(*pos)->set=true;
					pos++;
				} else {
					std::cerr << env.name << ":" << " invalid command -- '" << optopt << "' " << std::endl; 
					return -1; 
				}
			}
			
		};
	};
	if (env.commands.size()!=0) return -1;
	return env.required_set()!=1;
}

void 
print_version(Environment env){
	printf("%s version %s written by %s\n", env.name, env.version, env.author);
	exit(0);
}

void 
Usage(Environment env, FILE *out){
	std::string all_flags="";
	std::string required_args="";
	std::string optional_args="";
	size_t x=0;
	if (env.commands.size()!=0) fprintf(out, "usage: %s <command> [<args>]\n\n", env.name);
	else {
		std::list <Flag>::iterator flag=env.flags.begin();
		std::list <Flag>::iterator f_end=env.flags.end();
		while (flag!=f_end){
			all_flags+=flag->opt;
			flag++;
		}
		std::list <Argument>::iterator arg=env.args.begin();
		std::list <Argument>::iterator end=env.args.end();
		while (arg!=end){
			if (!arg->required) {
				if (x<6){
					optional_args+=" [--";
					optional_args+=arg->lopt;
					optional_args+=" ";
					optional_args+="VALUE";//arg->operand_type;
					optional_args+="]";
					if (x==3) {optional_args+=" ...";}
				}
			} else {
				required_args+=" <--";
				required_args+=arg->lopt;
				required_args+=" ";
				required_args+="VALUE";//arg->operand_type;
				required_args+=">";
			}
			++arg;
			++x;
		}
		optional_args+=required_args;
		if (env.flags.size()!=0) {
			fprintf(out, "usage: %s -%s %s\n\n", env.name, all_flags.c_str(), format_usage(optional_args.c_str(),8) );
		} else {
			fprintf(out, "usage: %s %s\n\n", env.name, format_usage(optional_args.c_str(),8) );
		}
	}
}

void 
print_help(Environment env)
{
	Usage(env, stdout);

	printf("%s version %s written by %s\n", env.name, env.version, env.author);
 	printf("%s\n\n", format_usage(env.description, 0) );

	std::list <Command>::iterator com=env.commands.begin();
	std::list <Command>::iterator com_end=env.commands.end();
	if (env.commands.size()!=0) printf("Commands:\n");
	while (com!=com_end){
		//TODO Better length formating.
		if (strlen(com->lopt)>5)
			printf("  %s\t\t%s\n", com->lopt, format_usage(com->umsg, 24) );
		else
			printf("  %s\t\t\t%s\n", com->lopt, format_usage(com->umsg, 24) );
		++com;
	}
	std::list <Argument>::iterator arg=env.args.begin();
	std::list <Argument>::iterator end=env.args.end();
	printf("Options:\n");
	while (arg!=end){
		//TODO Better length formating.
		if (strlen(arg->lopt)>7)
			printf("  -%c, --%s\t%s\n", arg->opt, arg->lopt, format_usage(arg->umsg, 24) );
		else
			printf("  -%c, --%s\t\t%s\n", arg->opt, arg->lopt, format_usage(arg->umsg, 24) );
		++arg;
	}
	std::list <Flag>::iterator flag=env.flags.begin();
	std::list <Flag>::iterator flag_end=env.flags.end();
	printf("Flags:\n");
	while (flag!=flag_end){
		//TODO Better length formating.
		if (strlen(flag->lopt)>7)
			printf("  -%c, --%s\t%s\n", flag->opt, flag->lopt, format_usage(flag->umsg, 24) );
		else
			printf("  -%c, --%s\t\t%s\n", flag->opt, flag->lopt, format_usage(flag->umsg, 24) );
		++flag;
	}
	printf("\n%s\n", format_usage(env.footer_, 0) );
	exit(0);
}


void 
print_usage(Environment env)
{
	Usage(env, stderr);
	fprintf(stderr, "Try '%s --help' for more information.\n", env.name);
	exit(0);
}

void 
Environment::close(void)
{
	flags.clear();
	args.clear();
	commands.clear();
	required_args.clear();
	required_coms.clear();
}

int
flag_commands(void *parm)
{
	print_commands(*(Environment *) parm);
	exit(0);
}

void print_commands(Environment env)
{
	std::list <Command>::iterator com=env.commands.begin();
	std::list <Command>::iterator com_end=env.commands.end();
	if (env.commands.size()==0) return;
	printf("%s", com->lopt);
	com++;
	while (com!=com_end){
		printf(" %s", com->lopt);
		com++;
	}
	printf("\n");
}

void
print_options(Environment env)
{
	std::list <Argument>::iterator arg=env.args.begin(), arg_end=env.args.end();
	std::list <Flag>::iterator flag=env.flags.begin(), flags_end=env.flags.end();

	if (env.args.size()!=0)
	{
		printf("%s", arg->lopt);
		arg++;
		while (arg!=arg_end) {printf(" %s", arg->lopt); arg++;}
	}
	if (env.flags.size()!=0)
	{
		printf("%s", flag->lopt);
		flag++;
		while (flag!=flags_end) {printf(" %s", flag->lopt); flag++;}
	}
	printf("\n");
	return;
}

int
flag_options(void *parm)
{
	print_options(*(Environment *) parm);
	exit(0);
}

void 
Environment::set_footer(const char *str)
{
	footer_=str;
}
