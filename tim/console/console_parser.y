%defines
%locations
%error-verbose
%name-prefix = "console_"

%{
// Description: GNU Bison-generated parser code
// Detail: Generated from console\_parser.y.
// Documentation: tim_console_cpp.txt

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

namespace Tim
{

	Program_AstNode* programAst;

}

namespace
{

	void setPosition(AstNode* node, const YYLTYPE& position)
	{
		node->setPosition(position.first_line, position.first_column);
	}

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

%token <string> T_INTEGER T_REAL T_IDENTIFIER T_PRINT T_STRING

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
		setPosition(programAst, @1);
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
		Print_AstNode* printNode = new Print_AstNode($2);
		setPosition(printNode, @1);
		$$ = printNode;
	}
	| expression
	{
		Declaration_AstNode* declaration = 
			new Declaration_AstNode("ans", $1);
		setPosition(declaration, @1);
		$$ = declaration;
	}
	;

declaration
	: identifier '=' expression
	{
		Declaration_AstNode* declaration = 
			new Declaration_AstNode(*$1, $3);
		setPosition(declaration, @1);
		$$ = declaration; 
	}
	;
	
function_call
	: identifier '(' expression_list ')'
	{
		FunctionCall_AstNode* functionCall = 
			new FunctionCall_AstNode(*$1, $3);
		setPosition(functionCall, @1);
		$$ = functionCall;
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
		Identifier_AstNode* identifier = 
			new Identifier_AstNode(*name);
		setPosition(identifier, @1);
		$$ = identifier;
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
		Integer_AstNode* integerNode = 
			new Integer_AstNode($1);
		setPosition(integerNode, @1);
		$$ = integerNode;
	}	
	| real_value
	{
		Real_AstNode* realNode = 
			new Real_AstNode($1);
		setPosition(realNode, @1);
		$$ = realNode;
	}
	| string
	{
		String_AstNode* stringNode =
			new String_AstNode(*$1);
		setPosition(stringNode, @1);
		$$ = stringNode;
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
		std::string* text = $1;
		const integer value = stringToInteger(*text);
		delete text;
		$$ = value;
	}
	;

real_value
	: T_REAL
	{
		std::string* text = $1;
		const real value = stringToReal(*text);
		delete text;
		$$ = value;
	}
	;

string
	: T_STRING
	{
		$$ = $1;
	}
	;

real_array
	: '[' real_array_content  ']'
	{
		RealSetSet* realSetSet = $2;
		
		// Find out the maximum width and height
		// of the matrix.

		const integer height = realSetSet->size();
		integer width = 0;
		for (integer y = 0;y < height;++y)
		{
			RealSet*& realSet = (*realSetSet)[y];
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
			RealSet* realSet = (*realSetSet)[y];
			std::copy(realSet->begin(), realSet->end(),
				signal->data().columnBegin(y));
			delete realSet;
		}
		delete realSetSet;
		
		RealArray_AstNode* realArray = 		
			new RealArray_AstNode(signal);
		setPosition(realArray, @1);
		$$ = realArray;
	}
	;

cell_array
	: '{' cell_array_content  '}'
	{
		StringSetSet* cellSetSet = $2;
		
		// Find out the maximum width and height
		// of the matrix.

		const integer height = cellSetSet->size();
		integer width = 0;
		for (integer y = 0;y < height;++y)
		{
			StringSet* cellSet = (*cellSetSet)[y];
			if (cellSet->size() > width)
			{
				width = cellSet->size();				
			}
		}
		
		// Copy the data to a string array.
		
		Array<std::string>* cellContent = new Array<std::string>(width, height);
		
		for (integer y = 0;y < height;++y)
		{
			StringSet* cellSet = (*cellSetSet)[y];
			std::copy(cellSet->begin(), cellSet->end(),
				cellContent->rowBegin(y));
			delete cellSet;
		}
		delete cellSetSet;
		
		CellArray_AstNode* cellArray = 		
			new CellArray_AstNode(cellContent);
		setPosition(cellArray, @1);
		$$ = cellArray;
	}
	;

real_array_content
	: real_list
	{
		RealSetSet* realSetSet = new RealSetSet;
		realSetSet->push_back($1);
		$$ = realSetSet;
	}
	| real_array_content ';' real_list
	{
		RealSetSet* realSetSet = $1;
		realSetSet->push_back($3);
		$$ = realSetSet;		
	}
	;

real_list
	: /* empty */
	{
		RealSet* realSet = new RealSet;
		$$ = realSet;
	}
	| real_list_1
	{
		$$ = $1;
	}
	;
			
real_list_1
	: number
	{
		RealSet* realSet = new RealSet;
		realSet->push_back($1);
		$$ = realSet;
	}
	| real_list_1 ',' number
	{
		RealSet* realSet = $1;
		realSet->push_back($3);
		$$ = realSet;
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
		$$ = stringSetSet;
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
		$$ = stringSet;
	}
	| cell_list_1 ',' identifier
	{
		std::string* name = $3;
		StringSet* stringSet = $1;
		stringSet->push_back(*name);
		delete name;
		$$ = stringSet;
	}
	;

number
	: integer_value
	{
		$$ = $1;
	}
	| real_value
	{
		$$ = $1;
	}
	;

%%

void console_error(const std::string& s)
{
	extern int console_lineno;
	extern char *console_text;
  
	std::cerr << "Line " << console_lineno << ": "
		<< "Unknown token '" << console_text << "'." << std::endl;
	if (!s.empty())
	{
		std::cerr << s << std::endl;
	}
}

void console_error(char *s)
{
	console_error(std::string(s));
}
