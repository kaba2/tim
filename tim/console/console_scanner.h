// Description: Tokenizer definitions
// Documentation: tim_console_cpp.txt

#ifndef TIM_CONSOLE_SCANNER_H
#define TIM_CONSOLE_SCANNER_H

#include "tim/console/ast.h"

#include <tim/core/mytypes.h>
#include <tim/core/signal.h>

#include <pastel/sys/array.h>

#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>

#include <vector>
#include <string>

namespace Tim
{

	typedef Array<SignalPtr> Cell;
	typedef boost::shared_ptr<Cell> CellPtr;

	typedef std::vector<std::string> StringSet;
	typedef std::vector<StringSet*> StringSetSet;

	typedef std::vector<real> RealSet;
	typedef std::vector<RealSet*> RealSetSet;

	typedef std::vector<boost::any> AnySet;

}

#include "tim/console/console_parser.tab.h"

#endif
