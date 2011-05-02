// Description: Main program
// Documentation: tim_console_cpp.txt

#include <iostream>
#include <string>

#include "tim/console/console_parser.h"
#include "tim/console/interpreter_astvisitor.h"
#include "tim/console/printer_astvisitor.h"
#include "tim/console/errorlog.h"

#include "pastel/sys/logging.h"

using namespace Tim;
using namespace Pastel;

int main(int argc, char **argv)
{
	log().addLogger(
		LoggerPtr(new Stream_Logger(&std::cout)));

	// Parse and construct an Abstract Syntax Tree.

	//std::cout << "Parsing..." << std::endl;
	//console_debug = 1;
	console_parse();

	// * If there is an error in the tokenizing phase,
	// the program will have already exited.
	// * If there is a syntax error in the parsing phase,
	// the program will continue here.

	if (programAst)
	{
		//Printer_AstVisitor printer(std::cout);
		//programAst->accept(printer);

		//std::cout << "Interpreting..." << std::endl;

		// Interpret the program.
		
		Interpreter_AstVisitor interpreter;
		try
		{
			programAst->accept(interpreter);
		}
		catch(const Interpreter_Exception&)
		{
			std::cerr << "Semantic errors found!" << std::endl;
		}

		delete programAst;
		programAst = 0;
	}

	// Output the possible error reports.

	std::cerr << errorLog();

	return 0;
}

