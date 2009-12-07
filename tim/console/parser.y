%defines
%locations
%error-verbose

%{
#define YY_NO_UNPUT
#define YYDEBUG 1

#include <iostream>
#include <string>
#include <stdio.h>
#include <map>

int yyerror(char *s);
int yylex(void);

#include <tim/core/mytypes.h>
#include <tim/core/signal_tools.h>

#include "tim/console/scanner.h"
#include "tim/console/errorlog.h"
#include "tim/console/functions.h"

#include <pastel/sys/stdext_copy_n.h>
#include <pastel/sys/string_tools.h>

#include <boost/any.hpp>
#include <boost/function.hpp>

using namespace Tim;

typedef std::map<std::string, boost::any*> SymbolMap;
typedef SymbolMap::const_iterator SymbolIterator;

struct YYLTYPE;

struct FunctionInfo
{
	typedef boost::function<boost::any*(const YYLTYPE& location, const AnySet& argSet)> Callback;

	FunctionInfo()
	: callback()
	, minArgs(0)
	, maxArgs(0)
	{
	}

	FunctionInfo(Callback callback_, integer minArgs_, integer maxArgs_)
	: callback(callback_)
	, minArgs(minArgs_)
	, maxArgs(maxArgs_)
	{
	}

	Callback callback;
	integer minArgs;
	integer maxArgs;
};

typedef std::map<std::string, FunctionInfo> FunctionMap;
typedef FunctionMap::const_iterator FunctionIterator;

bool initialized = false;

FunctionMap functionMap;
SymbolMap symbolMap;
ErrorLog errorLog;

namespace Tim
{

	void reportError(const YYLTYPE& location, const std::string& text);

	void print(boost::any* that);

	boost::any* functionCall(const YYLTYPE& location, 
		const std::string& name, const AnySet& argSet);

}

%}

%union{
	std::string* string;
	AnySet* anySet;
	RealSet* realSet;
	RealArray* realArray;
	CellSet* cellSet;
	CellArray* cellArray;
	Signal* signal;
	Cell* cell;
	real realValue;
	boost::any* any;
}

%token <string> T_INTEGER T_REAL T_IDENTIFIER T_PRINT T_GAUSSIAN

%type <realValue> number
%type <realSet> real_list real_list_1
%type <realArray> real_array_content
%type <signal> real_array signal_expression

%type <cellSet> cell_list cell_list_1
%type <cellArray> cell_array_content
%type <cell> cell_array 
%type <any> expression function_call
%type <anySet> expression_list expression_list_1

%left ';'
%left ','

%%

start
	: grammar
	;
	
grammar	
	: statement
	| grammar statement
	;

statement
	: declaration 
	| function_call
	| T_PRINT expression
	{
		print($2);
	}
	;

declaration
	: T_IDENTIFIER '=' expression
	{
		if (symbolMap.find(*$1) != symbolMap.end())
		{
			reportError(@1, "Multiple definitions for an identifier.");
		}
		symbolMap[*$1] = $3;
	}
	;
	
function_call
	: T_IDENTIFIER '(' expression_list ')'
	{
		$$ = functionCall(@1, *$1, *$3);
	}
	;
	
expression_list
	: /* empty */
	{
		$$ = new AnySet;
	}
	| expression_list_1
	{
		$$ = $1;
	}
	;
	
expression_list_1
	: expression
	{
		$$ = new AnySet;
		$$->push_back($1);
	}
	| expression_list_1 ',' expression
	{
		$$ = $1;
		$$->push_back($3);
	}
	;

expression
	: T_IDENTIFIER
	{
		SymbolIterator iter = symbolMap.find(*$1);
		if (iter == symbolMap.end())
		{
			$$ = new boost::any;
			reportError(@1, "Undefined identifier.");
		}
		else
		{
			$$ = iter->second;
		}
	}
	| function_call
	{
		$$ = $1;
	}
	| real_array
	{
		$$ = new boost::any($1);
	}
	| cell_array
	{
		$$ = new boost::any($1);
	}
	| T_INTEGER
	{
		$$ = new boost::any(stringToInteger(*$1));
	}	
	| T_REAL
	{
		$$ = new boost::any(stringToReal(*$1));
	}	
	;

signal_expression
	: T_IDENTIFIER
	{
		SymbolIterator iter = symbolMap.find(*$1);
		if (iter == symbolMap.end())
		{
			$$ = (Signal*)0;
			reportError(@1, "Undefined identifier.");
		}
		else
		{
			try
			{
				$$ = boost::any_cast<Signal*>(*iter->second);
			}
			catch(const boost::bad_any_cast&)
			{
				reportError(@1, "Identifier is not a real array.");
			}
		}
	}
	| real_array
	{
		$$ = $1;
	}
	;

real_array
	: T_GAUSSIAN
	{
		SignalPtr signal = generateGaussian(10000, 10);
		$$ = new Signal(*signal);
	}
	| '[' real_array_content  ']'
	{
		RealArray* realArray = $2;
		
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

		$$ = new Signal(samples, dimension);
		
		for (integer y = 0;y < height;++y)
		{
			RealSet*& realSet = (*realArray)[y];
		
			StdExt::copy_n(
				realSet->begin(),
				std::min((integer)realSet->size(), samples),
				$$->data().columnBegin(y));
			
			delete realSet;
			realSet = 0;
		}
		delete realArray;
	}
	;

cell_array
	: '{' cell_array_content  '}'
	{
		CellArray* cellArray = $2;
		
		// Find out the maximum width and height
		// of the matrix.

		const integer height = cellArray->size();
		integer width = 0;
		for (integer y = 0;y < height;++y)
		{
			CellSet*& cellSet = (*cellArray)[y];
			if (cellSet->size() > width)
			{
				width = cellSet->size();				
			}
		}
		
		// Copy the data to a cell array.
		
		$$ = new Cell(width, height);
		
		for (integer y = 0;y < height;++y)
		{
			CellSet*& cellSet = (*cellArray)[y];
			
			const integer cells = std::min((integer)cellSet->size(), width);
			for (integer x = 0;x < cells;++x)
			{
				(*$$)(x, y) = SignalPtr((*cellSet)[x]);
			}
			
			delete cellSet;
			cellSet = 0;
		}
		delete cellArray;
	}
	;

real_array_content
	: real_list
	{
		$$ = new RealArray;
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
		$$ = new CellArray;
		$$->push_back($1);
	}
	| cell_array_content ';' cell_list
	{
		$$ = $1;
		$$->push_back($3);
	}
	;

cell_list
	: /* empty */
	{
		$$ = new CellSet;
	}
	| cell_list_1
	{
		$$ = $1;
	}
	;
	
cell_list_1
	: signal_expression
	{
		$$ = new CellSet;
		$$->push_back($1);
	}
	| cell_list_1 ',' signal_expression
	{
		$$ = $1;
		$$->push_back($3);
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

int yyerror(const std::string& s)
{
	extern int yylineno;
	extern char *yytext;
  
	std::cerr << "ERROR: " << s << " at symbol \"" << yytext;
	std::cerr << "\" on line " << yylineno << std::endl;
	exit(1);
}

int yyerror(char *s)
{
	return yyerror(std::string(s));
}

namespace Tim
{

	void reportError(const YYLTYPE& location, const std::string& text)
	{
		errorLog.report(location.last_line, text);
	}
	
	void printErrors(std::ostream& stream)
	{
		stream << errorLog;	
	}

	void print(boost::any* that)
	{
		try
		{
			Signal* signal = boost::any_cast<Signal*>(*that);
			std::cout << *signal << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			integer k = boost::any_cast<integer>(*that);
			std::cout << k << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}

		try
		{
			real k = boost::any_cast<real>(*that);
			std::cout << k << std::endl;
		}
		catch(const boost::bad_any_cast&)
		{
		}
	}


	boost::any* functionCall(const YYLTYPE& location, 
		const std::string& name, const AnySet& argSet)
	{
		if (!initialized)
		{
			functionMap.insert(
				std::make_pair("differential_entropy_kl", 
				FunctionInfo(differential_entropy_kl, 1, 3)));
			functionMap.insert(
				std::make_pair("differential_entropy_kl_t", 
				FunctionInfo(differential_entropy_kl_t, 2, 4)));
				
			functionMap.insert(
				std::make_pair("differential_entropy_nk", 
				FunctionInfo(differential_entropy_nk, 1, 2)));

			functionMap.insert(
				std::make_pair("divergence_wkv", 
				FunctionInfo(divergence_wkv, 2, 2)));

			functionMap.insert(
				std::make_pair("mutual_information_t", 
				FunctionInfo(mutual_information_t, 3, 6)));
			functionMap.insert(
				std::make_pair("mutual_information", 
				FunctionInfo(mutual_information, 2, 5)));
			functionMap.insert(
				std::make_pair("mutual_information_pt", 
				FunctionInfo(mutual_information_pt, 4, 8)));
			functionMap.insert(
				std::make_pair("mutual_information_p", 
				FunctionInfo(mutual_information_p, 3, 7)));
				
			functionMap.insert(
				std::make_pair("transfer_entropy_t", 
				FunctionInfo(transfer_entropy_t, 4, 8)));
			functionMap.insert(
				std::make_pair("transfer_entropy_pt", 
				FunctionInfo(transfer_entropy_pt, 5, 10)));
			functionMap.insert(
				std::make_pair("transfer_entropy", 
				FunctionInfo(transfer_entropy, 3, 7)));
			functionMap.insert(
				std::make_pair("transfer_entropy_p", 
				FunctionInfo(transfer_entropy_p, 4, 9)));

			initialized = true;
		}
	
		FunctionIterator iter = functionMap.find(name);
		if (iter == functionMap.end())
		{
			reportError(location, "Undefined function.");
		}
		else
		{		
			const FunctionInfo& info = iter->second;
			bool error = false;
			if (argSet.size() > info.maxArgs)
			{
				reportError(location, "Too many arguments to a function (max " + integerToString(info.maxArgs) + ").");
				error = true;
			}
			if (argSet.size() < info.minArgs)
			{
				reportError(location, "Not enough arguments to a function (min " + integerToString(info.minArgs) + ").");
				error = true;
			}
			
			if (!error)
			{
				return info.callback(location, argSet);
			}
		}
		
		return 0;
	}

}
