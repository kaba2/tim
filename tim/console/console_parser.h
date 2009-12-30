// Description: Parser and interpreter
// Documentation: tim_console_cpp.txt

#ifndef TIM_CONSOLE_PARSER_H
#define TIM_CONSOLE_PARSER_H

#include "tim/console/console_scanner.h"
#include "tim/console/ast.h"

#include <iostream>
#include <string>

int console_parse();
extern int console_debug;

namespace Tim
{

	void printErrors(std::ostream& stream);
	void reportError(const YYLTYPE& location, const std::string& text);

	extern Program_AstNode* programAst;

}

#endif
