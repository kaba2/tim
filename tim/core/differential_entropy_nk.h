// Description: Differential entropy estimation
// Detail: Nilsson-Kleijn manifold nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_H
#define TIM_DIFFERENTIAL_ENTROPY_NK_H

#include "tim/core/mytypes.h"

#include <pastel/sys/range.h>

namespace Tim
{

	template <
		typename Signal_Iterator, 
		typename NormBijection>
	real differentialEntropyNk(
		const boost::iterator_range<Signal_Iterator>& signalSet,
		const NormBijection& normBijection,
		integer* outIntrinsicDimension = 0);

}

#include "tim/core/differential_entropy_nk.hpp"

#endif
