// Description: Type definitions for the TIM library
// Detail: Simply imports the type definitions of the Pastel library
// Documentation: basics.txt

#ifndef TIM_MYTYPES_H
#define TIM_MYTYPES_H

#include "tim/core/tim_library.h"

#include <pastel/sys/mytypes.h>
#include <pastel/sys/vector_tools.h>

#include <pastel/math/maximum_normbijection.h>
#include <pastel/geometry/splitrules.h>

#if (defined _WIN32 || defined _WIN64)
#	pragma detect_mismatch("TIM_VERSION", "1.3")
#endif

namespace Tim
{

	using namespace Pastel;

	typedef Maximum_NormBijection<real> Default_NormBijection;
	typedef SlidingMidpoint_SplitRule SplitRule;

}

#endif
