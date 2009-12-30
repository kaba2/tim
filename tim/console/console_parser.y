%defines
%locations
%error-verbose
%name-prefix = "console_"

%{
#include <iostream>
#include <string>
#include <stdio.h>
#include <map>

#include <tim/core/mytypes.h>
#include <tim/core/signal_tools.h>

#include "tim/console/console_scanner.h"
#include "tim/console/errorlog.h"
#include "tim/console/functions.h"

#include <pastel/sys/stdext_copy_n.h>
#include <pastel/sys/string_tools.h>

#include <boost/any.hpp>
#include <boost/function.hpp>

using namespace Tim;

void console_error(char *s);
int console_lex();

ErrorLog errorLog;

namespace Tim
{

	void reportError(const YYLTYPE& location, const std::string& text);

	Program_AstNode* programAst;

}

%}

%union{
	Tim::Statement_AstNode* statement;
	Tim::Expression_AstNode* expression;
	Tim::Declaration_AstNode* declaration;
	Tim::FunctionCall_AstNode* function_call;
	Tim::CellArray_AstNode* cell_array;
	Tim::RealArray_AstNode* real_array;
	Tim::StatementSet* program;
	
	Tim::StringSet* cell_list;
	Tim::StringSetSet* cell_array_content;
	Tim::RealSet* real_list;
	Tim::RealSetSet* real_array_content;
	Tim::ExpressionSet* expression_list;
	Tim::real real_value;
	Tim::integer integer_value;

	std::string* string;
}

%token <string> T_INTEGER T_REAL T_IDENTIFIER T_PRINT T_GAUSSIAN T_STRING

%type <program> program
%type <statement> statement
%type <expression> expression
%type <declaration> declaration
%type <function_call> function_call
%type <cell_array> cell_array
%type <cell_list> cell_list cell_list_1
%type <cell_array_content> cell_array_content

%type <real_array> real_array
%type <real_array_content> real_array_content
%type <real_list> real_list real_list_1

%type <expression_list> expression_list expression_list_1

%type <string> identifier string
%type <integer_value> integer_value
%type <real_value> number real_value

%%

start
	: program
	{
		programAst = new Program_AstNode($1);		
	}
	;
	
program
	: statement
	{
		StatementSet* statementSet = new StatementSet;
		statementSet->push_back($1);
		$$ = statementSet;
	}
	| program statement
	{
		StatementSet* statementSet = $1;
		statementSet->push_back($2);
		$$ = statementSet;
	}
	;

statement
	: declaration 
	{
		$$ = $1;
	}
	| T_PRINT expression
	{
		$$ = new Print_AstNode($2);
	}
	| expression
	{
		$$ = new Declaration_AstNode("ans", $1);
	}
	;

declaration
	: identifier '=' expression
	{
		$$ = new Declaration_AstNode(*$1, $3);
	}
	;
	
function_call
	: identifier '(' expression_list ')'
	{
		$$ = new FunctionCall_AstNode(*$1, $3);
	}
	;
	
expression_list
	: /* empty */
	{
		ExpressionSet* expressionSet = new ExpressionSet;
		$$ = expressionSet;
	}
	| expression_list_1
	{
		ExpressionSet* expressionSet = $1;
		$$ = expressionSet;
	}
	;
	
expression_list_1
	: expression
	{
		ExpressionSet* expressionSet = new ExpressionSet;
		expressionSet->push_back($1);
		$$ = expressionSet;
	}
	| expression_list_1 ',' expression
	{
		ExpressionSet* expressionSet = $1;
		expressionSet->push_back($3);
		$$ = expressionSet;
	}
	;

expression
	: identifier
	{
		std::string* name = $1;
		$$ = new Identifier_AstNode(*$1);
		delete name;
	}
	| function_call
	{
		$$ = $1;
	}
	| real_array
	{
		$$ = $1;
	}
	| cell_array
	{
		$$ = $1;
	}
	| integer_value
	{
		$$ = new Integer_AstNode($1);
	}	
	| real_value
	{
		$$ = new Real_AstNode($1);
	}
	| string
	{
		$$ = new String_AstNode(*$1);
	}
	;

identifier
	: T_IDENTIFIER
	{
		$$ = $1;
	}
	;
	
integer_value
	: T_INTEGER
	{
		$$ = stringToInteger(*$1);
	}
	;

real_value
	: T_REAL
	{
		$$ = stringToReal(*$1);
	}
	;

string
	: T_STRING
	{
		$$ = $1;
	}
	;

real_array
	: T_GAUSSIAN
	{
		$$ = new RealArray_AstNode(generateGaussian(10000, 10));
	}
	| '[' real_array_content  ']'
	{
		RealSetSet* realArray = $2;
		
		// Find out the maximum width and height
		// of the matrix.

		const integer height = realArray->size();
		integer width = 0;
		for (integer y = 0;y < height;++y)
		{
			RealSet*& realSet = (*realArray)[y];
			if (realSet->size() > width)
			{
				width = realSet->size();				
			}
		}
		
		// Copy the data to a signal.
		
		const integer samples = width;
		const integer dimension = height;

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));
		
		for (integer y = 0;y < height;++y)
		{
			RealSet* realSet = (*realArray)[y];
			std::copy(realSet->begin(), realSet->end(),
				signal->data().columnBegin(y));
			delete realSet;
		}
		delete realArray;
		
		$$ = new RealArray_AstNode(signal);
	}
	;

cell_array
	: '{' cell_array_content  '}'
	{
		StringSetSet* cellArray = $2;
		
		// Find out the maximum width and height
		// of the matrix.

		const integer height = cellArray->size();
		integer width = 0;
		for (integer y = 0;y < height;++y)
		{
			StringSet*& cellSet = (*cellArray)[y];
			if (cellSet->size() > width)
			{
				width = cellSet->size();				
			}
		}
		
		// Copy the data to a string array.
		
		Array<std::string>* cellContent = new Array<std::string>(width, height);
		
		for (integer y = 0;y < height;++y)
		{
			StringSet* cellSet = (*cellArray)[y];
			std::copy(cellSet->begin(), cellSet->end(),
				cellContent->rowBegin(y));
			delete cellSet;
		}
		delete cellArray;
		
		$$ = new CellArray_AstNode(cellContent);
	}
	;

real_array_content
	: real_list
	{
		$$ = new RealSetSet;
		$$->push_back($1);
	}
	| real_array_content ';' real_list
	{
		$$ = $1;
		$$->push_back($3);
	}
	;

real_list
	: /* empty */
	{
		$$ = new RealSet;
	}
	| real_list_1
	{
		$$ = $1;
	}
	;
			
real_list_1
	: number
	{
		$$ = new RealSet;
		$$->push_back($1);
	}
	| real_list_1 ',' number
	{
		$$ = $1;
		$$->push_back($3);
	}
	;

cell_array_content
	: cell_list
	{
		StringSetSet* stringSetSet = new StringSetSet;
		stringSetSet->push_back($1);
		$$ = stringSetSet;
	}
	| cell_array_content ';' cell_list
	{
		StringSetSet* stringSetSet = $1;
		stringSetSet->push_back($3);
	}
	;

cell_list
	: /* empty */
	{
		StringSet* stringSet = new StringSet;
		$$ = stringSet;
	}
	| cell_list_1
	{
		StringSet* stringSet = $1;
		$$ = stringSet;
	}
	;
	
cell_list_1
	: identifier
	{
		std::string* name = $1;
		StringSet* stringSet = new StringSet;
		stringSet->push_back(*name);
		delete name;
	}
	| cell_list_1 ',' identifier
	{
		std::string* name = $3;
		StringSet* stringSet = $1;
		stringSet->push_back(*name);
		delete name;
	}
	;

number
	: T_INTEGER
	{
		$$ = stringToInteger(*$1);
	}
	| T_REAL
	{
		$$ = stringToReal(*$1);
	}
	;

%%

void console_error(const std::string& s)
{
	extern int console_lineno;
	extern char *console_text;
  
	std::cerr << "ERROR: " << s << " at symbol \"" << console_text;
	std::cerr << "\" on line " << console_lineno << std::endl;
}

void console_error(char *s)
{
	console_error(std::string(s));
}

namespace Tim
{

	void reportError(const YYLTYPE& location, const std::string& text)
	{
		errorLog.report(location.first_line, text);
	}
	
	void printErrors(std::ostream& stream)
	{
		stream << errorLog;	
	}

}
