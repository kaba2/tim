%{
#include <iostream>
#include <stdio.h>
#include <string>

#include "tim/console/scanner.h"

using namespace Tim;

#include "parser.tab.h"
void yyerror(char *s);
int yylineno = 1;
%}

%option never-interactive nounput

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
				return yytext[0]; 
			}

{IDENTIFIER}	{
					yylval.string = new std::string(yytext);
					return T_IDENTIFIER;
				}

{INTEGER}	{
				yylval.string = new std::string(yytext);
				return T_INTEGER; 
			}

{REAL}		{ 
				yylval.string = new std::string(yytext);
				return T_REAL; 
			}

{STRING}	{
				std::string copy(yytext);
				yylval.string = new std::string(copy.substr(1, copy.size() - 2));
				return T_STRING;
			}

[ \t\r]*	{
				/* Eat whitespace (new-lines are handled specially). */
			}

[\n]		{
				/* Eat new-lines but count rows. */
				++yylineno;
			}

.			{ 
				std::cerr << "SCANNER "; 
				yyerror(""); 
				exit(1);	
			}


