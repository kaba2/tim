// Description: Main program
// Documentation: tim_console_cpp.txt

/*
Syntax
------

A context-free grammar for the TimConsole's script files. 

As an abbreviation we use:
non-terminal % separator = (non-terminal (separator non-terminal)*)?

grammar := statement*
statement := declaration | function_call | 'print' expression
declaration := identifier '=' expression 
expression := identifier | real_array | cell_array | integer | real
function_call := identifier '(' expression % ',' ')'

signal_expression := identifier | real_array
real_array := '[' real_list % ';' ']'
real_list := real % ','

cell_expression := identifier | cell_array
cell_array := '{' cell_list % ';' '}'
cell_list := signal_expression % ','

identifier := alpha (alpha | digit)*

Syntax example
--------------

integer k = 1
real maxRelativeError = 2

RealArray A = [3 x 3][1.1111, 2,     3;
4,      +0.05, 6;
7,      -8,    9]
RealArray B = A
CellArray C = {2 x 2}{A, A; B, B}

RealArray S = differentialEntropyKl({A})
print(S)

Semantics
---------

- An identifier must be declared before use.
- Row-comments are given by // or %.
- Comments are ignored.
- Tokens are separated by white-space (space, tab, newline). 
- The amount and type of white-space is irrelevant. I.e. free form.
*/

#include <iostream>
#include <string>

#include "tim/console/parser.h"

using namespace Tim;

int main(int argc, char **argv)
{
	//yydebug = 1;
	yyparse();

	printErrors(std::cerr);

	return 0;
}

