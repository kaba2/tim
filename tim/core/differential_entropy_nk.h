// Description: Differential entropy estimation
// Detail: Nilsson-Kleijn manifold nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_H
#define TIM_DIFFERENTIAL_ENTROPY_NK_H

#include "tim/core/mytypes.h"

#include <pastel/sys/range.h>

namespace Tim
{

	template <
		typename SignalPtr_Range, 
		typename NormBijection>
	real differentialEntropyNk(
		const SignalPtr_Range& signalSet,
		const NormBijection& normBijection,
		integer* outIntrinsicDimension = 0);

}

#include "tim/core/differential_entropy_nk.hpp"

#endif
