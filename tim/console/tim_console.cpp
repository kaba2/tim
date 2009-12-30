// Description: Main program
// Documentation: tim_console_cpp.txt

#include <iostream>
#include <string>

#include "tim/console/console_parser.h"
#include "tim/console/interpreter.h"

using namespace Tim;

int main(int argc, char **argv)
{
	// Parse and construct an Abstract Syntax Tree.

	//console_debug = 1;
	console_parse();

	// Interpret the program.
	
	Interpreter_AstVisitor interpreter;
	programAst->accept(interpreter);

	// Output the possible error reports.

	printErrors(std::cerr);

	return 0;
}

