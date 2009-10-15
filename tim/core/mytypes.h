// Description: Type definitions for the Tim library
// Detail: Simply imports the type definitions of the Pastel library
// Documentation: basics.txt

#ifndef TIM_MYTYPES_H
#define TIM_MYTYPES_H

#include "tim/core/tim_library.h"

#include <pastel/sys/mytypes.h>
#include <pastel/sys/vector_tools.h>

#include <pastel/math/maximum_normbijection.h>

namespace Tim
{

	using namespace Pastel;

	typedef Maximum_NormBijection<real> Default_NormBijection;

}

#endif
