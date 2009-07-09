#ifndef TIM_MUTUAL_INFORMATION_NAIVE_H
#define TIM_MUTUAL_INFORMATION_NAIVE_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

#include <pastel/sys/smallset.h>

namespace Tim
{

	TIMCORE void mutualInformationNaive(
		const SignalPtr& signal,
		integer bins,
		MatrixD& result);

	template <typename NormBijection>
	real mutualInformationFromEntropy(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

}

#endif
