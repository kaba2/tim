// Description: Differential entropy estimation
// Detail: Nilsson-Kleijn manifold nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_H
#define TIM_DIFFERENTIAL_ENTROPY_NK_H

#include "tim/core/mytypes.h"

#include <pastel/sys/iteratorrange.h>

namespace Tim
{

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropyNk(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		real maxRelativeError,
		const NormBijection& normBijection,
		integer* outIntrinsicDimension = 0);

}

#include "tim/core/differential_entropy_nk.hpp"

#endif
