%option never-interactive nounput noyywrap
%option prefix = "console_"
%option yylineno

%{
#include <iostream>
#include <string>

#define YY_USER_ACTION \
	{ \
		console_lloc.first_line = yylineno; \
		console_lloc.last_line = yylineno; \
	}

#include "tim/console/console_scanner.h"

void console_error(char *s);
%}

ALPHA		[A-Z]|[a-z]|_
DIGIT		[0-9]
ALNUM		{ALPHA}|{DIGIT}
IDENTIFIER	{ALPHA}{ALNUM}*
INTEGER		[+-]?[0-9]+
REAL		[+-]?[0-9]*\.[0-9]+([Ee][+-]?[0-9]+)?
SHORTIES	[,;\[\]\*={}()]
STRING		\"[^\"]*\"

%%

\/\/.*      {
				/* Eat C-style one-line comments */
			}

%.*         {
				/* Eat Matlab-style one-line comments */
			}

print		{
				return T_PRINT;
			}

gaussian	{
				return T_GAUSSIAN;
			}

{SHORTIES}	{ 
				return console_text[0]; 
			}

{IDENTIFIER}	{
					console_lval.string = new std::string(console_text);
					return T_IDENTIFIER;
				}

{INTEGER}	{
				console_lval.string = new std::string(console_text);
				return T_INTEGER; 
			}

{REAL}		{ 
				console_lval.string = new std::string(console_text);
				return T_REAL; 
			}

{STRING}	{
				std::string copy(console_text);
				console_lval.string = new std::string(copy.substr(1, copy.size() - 2));
				return T_STRING;
			}

[ \t\r\n]*	{
				/* Eat whitespace. */
			}

.			{ 
				console_error(""); 
				exit(1);	
			}


