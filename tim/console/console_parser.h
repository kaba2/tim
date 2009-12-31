// Description: Parser definitions
// Detail: GNU Bison is used to generate a scanner from a grammar.
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

	extern Program_AstNode* programAst;

}

#endif
