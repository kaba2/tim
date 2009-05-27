#ifndef TIM_DIFFERENTIAL_ENTROPY_H
#define TIM_DIFFERENTIAL_ENTROPY_H

#include "tim/core/signal.h"

namespace Tim
{

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

}

#include "tim/core/differential_entropy.hpp"

#endif
