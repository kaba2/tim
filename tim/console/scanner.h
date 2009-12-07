#ifndef TIM_SCANNER_H
#define TIM_SCANNER_H

#include <tim/core/mytypes.h>
#include <tim/core/signal.h>

#include <pastel/sys/array.h>

#include <boost/any.hpp>

#include <vector>

namespace Tim
{

	typedef Array<SignalPtr> Cell;

	typedef std::vector<real> RealSet;
	typedef std::vector<RealSet*> RealArray;

	typedef std::vector<Signal*> CellSet;
	typedef std::vector<CellSet*> CellArray;

	typedef std::vector<boost::any*> AnySet;

}

using namespace Tim;

#include "tim/console/parser.tab.h"

#endif
