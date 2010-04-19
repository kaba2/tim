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
INTEGER		[+-]?[0-9]+
REAL_NAN	[Nn][Aa][Nn]
REAL_INF	[Ii][Nn][Ff]
REAL_NUMBER [0-9]*\.[0-9]+([Ee][+-]?[0-9]+)?
REAL		[+-]?({REAL_NAN}|{REAL_INF}|{REAL_NUMBER})
SHORTIES	[,;\[\]\*={}()]
STRING		\"[^\"]*\"
IDENTIFIER	{ALPHA}{ALNUM}*

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

{SHORTIES}	{ 
				return console_text[0]; 
			}

{INTEGER}	{
				console_lval.string = new std::string(console_text);
				return T_INTEGER; 
			}

{REAL}		{ 
				console_lval.string = new std::string(console_text);
				return T_REAL; 
			}

{IDENTIFIER}	{
					console_lval.string = new std::string(console_text);
					return T_IDENTIFIER;
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


