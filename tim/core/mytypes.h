// Description: Type definitions for the Tim library
// Detail: Simply imports the type definitions of the Pastel library
// Documentation: basics.txt

#ifndef TIM_MYTYPES_H
#define TIM_MYTYPES_H

#include "tim/core/tim_library.h"

#include <pastel/sys/mytypes.h>
#include <pastel/sys/vector_tools.h>

#include <pastel/math/infinity_normbijection.h>

namespace Tim
{

	using namespace Pastel;

	typedef Infinity_NormBijection<real> Default_NormBijection;

}

#endif
