// Description: Parser and interpreter
// Documentation: tim_console_cpp.txt

#ifndef TIM_PARSER_H
#define TIM_PARSER_H

#include "tim/console/scanner.h"

#include <iostream>
#include <string>

int yyparse();
extern int yydebug;

namespace Tim
{

	void printErrors(std::ostream& stream);
	void reportError(const YYLTYPE& location, const std::string& text);

}

#endif
