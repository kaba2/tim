// Description: Built-in functions for the interpreter.

#ifndef TIM_FUNCTIONS_H
#define TIM_FUNCTIONS_H

#include "tim/console/console_scanner.h"

#include <string>

#include <boost/any.hpp>

namespace Tim
{

	class FunctionCall_Exception
	{
	};

	boost::any functionCall(const std::string& name, const AnySet& argSet);

}

#endif
