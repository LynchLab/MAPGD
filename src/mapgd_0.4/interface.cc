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
	//	*word='\n';
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
arg_setvectorstr(int argc, char **argv, void *parm)
{
	int n=1;
	std::vector <std::string> *v=(std::vector <std::string> *)(parm);
	while (n<argc){
	//	std::cout << "n\n";
		if (argv[n][0]=='-') return n;
	//	std::cout << "pusshing " << argv[n] << std::endl;
		v->push_back(argv[n]);
		if ( n+1==argc) return n+1;
		++n;
	} 
	std::cerr << __FILE__ << ":" << __LINE__ << ": error parsing " << argv[1] << std::endl;
	exit(1);
}

/*@Breif : sets a vector of ints from a string. */
int arg_setvectorint(int argc, char **argv, void *parm)
{
	std::vector <size_t> *v=(std::vector <size_t> *)(parm);
	if (argc>1){
		std::vector<std::string> elems=split(argv[1], ',');
		for (size_t x=0; x<elems.size(); ++x){
			if ( isint(elems[x].c_str() ) ) v->push_back(atoi(elems[x].c_str() )-1) ;
			else {
				std::vector<std::string> intpair=split(elems[x], '-');
				if (intpair.size()==2){
					if (isint(intpair[0].c_str() ) && isint(intpair[1].c_str() ) ){
						for (int y=atoi(intpair[0].c_str() ); y<=atoi(intpair[1].c_str() ); ++y){
							if (std::find(v->begin(), v->end(), y-1)==v->end() ){
								v->push_back(y-1);
							}
						};
					} else {
						std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
					};
				} else{
					std::cerr << "cannot parse string " << elems[x] << " into exactly two integers. Please use x-y formating." << std::endl;
				};
			}
		}
		return 2;
	} 
	std::cerr << "arg_setvectorint: error parsing " << argv[1] << std::endl;
	exit(1);
}

//TODO Comment
int arg_setint(int argc, char **argv, void *parm)
{
	if (argc>1){
		if (isint(argv[1])){
			*( (int *)(parm) )=atoi(argv[1]);
			return 2;
		} return argc+1;
	} 
	std::cerr << "arg_setint:error parsing " << argv[1] << std::endl;
	exit(1);
}

//TODO Comment
int arg_set2int(int argc, char **argv, void *parm)
{
	if (argc>2){
		if (isint(argv[1]) && isint(argv[2]) ){
			( (int *)(parm) )[0]=atoi(argv[1]);
			( (int *)(parm) )[1]=atoi(argv[2]);
			return 3;
		} return argc+1;
	} exit(1);
}


//TODO Comment
int arg_setchar(int argc, char **argv, void *parm)
{
	if (argc>1){
		*(char *)(parm)=argv[1][0];
		return 2;
	} exit(1);	
}

//TODO Comment
int arg_setstr(int argc, char **argv, void *parm)
{
        if (argc>1){
                *(std::string *)(parm)=argv[1];
                return 2;
        } exit(1);
}

//TODO Comment
int arg_setc_str(int argc, char **argv, void *parm)
{
	if (argc>1){
		*(char **)(parm)=argv[1];
		return 2;
	} exit(1);	
}

//TODO Comment
int arg_setfloat_t(int argc, char **argv, void *parm)
{
	if (argc>1){
		*(float_t *)(parm)=atof(argv[1]);
		return 2;
	} exit(1);	
}


//TODO Comment
int flag_set(void *parm)
{
	*(bool *)(parm)=!(*(bool *)(parm));
	return 0;
}

//TODO Comment
int flag_version(void *parm)
{
	print_version(*(env_t *) parm);
	return 0;
}


//! prints the long help statment to stdout (print_help).
int flag_help(void *parm)
{
	print_help(*(env_t *) parm);
	return 0;
}

//TODO Comment
int flag_usage(void *parm)
{
	print_usage(*(env_t *) parm);
	return 0;
}

//TODO Comment
int arg_error(int argc, char **argv, void *parm)
{
	std::cerr << "error: option functor unset\n";
	return 1;
}

/* The main command line argument parsing routine */

int 
parsargs(int argc, char *argv[], env_t &env)
{
	int optind=1; 
	char *optopt; 

	std::list <arg_t>::iterator arg=env.args.begin();
	std::list <arg_t>::iterator arg_end=env.args.end();

	std::list <flag_t>::iterator flag=env.flags.begin();
	std::list <flag_t>::iterator flag_end=env.flags.end();

	std::list <com_t>::iterator com=env.commands.begin();
	std::list <com_t>::iterator com_end=env.commands.end();

	while (optind<argc){
		if (argv[optind][0]=='-'){
			if (argv[optind][1]!=0){
				if (argv[optind][1]!='-'){
					optopt=argv[optind]+1;
					arg=env.args.begin();
					while(arg!=arg_end){
						if (*optopt==arg->opt){
							optind+=arg->func(argc-optind, argv+optind, arg->parm);
							arg->set=true;
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
				std::cerr << env.name << ":" << " invalid command -- '" << optopt << "' " << std::endl; 
				return -1; 
			}
			
		};
	};
	if (env.commands.size()!=0) return -1;
	return env.required_set()!=1;
}

void print_version(env_t env){
	printf("%s version %s written by %s\n", env.name, env.version, env.author);
	exit(0);
}

void Usage(env_t env, FILE *out){
	std::string all_flags="";
	std::string required_args="";
	std::string optional_args="";
	size_t x=0;
	if (env.commands.size()!=0) fprintf(out, "usage: %s <command> [<args>]\n\n", env.name);
	else {
		std::list <flag_t>::iterator flag=env.flags.begin();
		std::list <flag_t>::iterator f_end=env.flags.end();
		while (flag!=f_end){
			all_flags+=flag->opt;
			flag++;
		}
		std::list <arg_t>::iterator arg=env.args.begin();
		std::list <arg_t>::iterator end=env.args.end();
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
			fprintf(out, "usage: %s %s\n\n", env.name, all_flags.c_str(), format_usage(optional_args.c_str(),8) );
		}
	}
}

void 
print_help(env_t env)
{
	Usage(env, stdout);

	printf("%s version %s written by %s\n", env.name, env.version, env.author);
 	printf("%s\n\n", format_usage(env.description, 0) );

	std::list <com_t>::iterator com=env.commands.begin();
	std::list <com_t>::iterator com_end=env.commands.end();
	if (env.commands.size()!=0) printf("Commands:\n");
	while (com!=com_end){
		//TODO Better length formating.
		if (strlen(com->lopt)>5)
			printf("  %s\t\t%s\n", com->lopt, format_usage(com->umsg, 24) );
		else
			printf("  %s\t\t\t%s\n", com->lopt, format_usage(com->umsg, 24) );
		++com;
	}
	std::list <arg_t>::iterator arg=env.args.begin();
	std::list <arg_t>::iterator end=env.args.end();
	printf("Options:\n");
	while (arg!=end){
		//TODO Better length formating.
		if (strlen(arg->lopt)>7)
			printf("  -%c, --%s\t%s\n", arg->opt, arg->lopt, format_usage(arg->umsg, 24) );
		else
			printf("  -%c, --%s\t\t%s\n", arg->opt, arg->lopt, format_usage(arg->umsg, 24) );
		++arg;
	}
	std::list <flag_t>::iterator flag=env.flags.begin();
	std::list <flag_t>::iterator flag_end=env.flags.end();
	printf("Flags:\n");
	while (flag!=flag_end){
		//TODO Better length formating.
		if (strlen(flag->lopt)>7)
			printf("  -%c, --%s\t%s\n", flag->opt, flag->lopt, format_usage(flag->umsg, 24) );
		else
			printf("  -%c, --%s\t\t%s\n", flag->opt, flag->lopt, format_usage(flag->umsg, 24) );
		++flag;
	}
	exit(0);
}


void print_usage(env_t env){
	Usage(env, stderr);
	fprintf(stderr, "Try '%s --help' for more information.\n", env.name);
	exit(0);
}

void env_t::close(void){
	flags.clear();
	args.clear();
	commands.clear();
	required_args.clear();
	required_coms.clear();
}
