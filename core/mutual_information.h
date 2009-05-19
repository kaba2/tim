#ifndef TIM_MUTUAL_INFORMATION_H
#define TIM_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

namespace Tim
{

	//! Computes mutual information between a set of signals.

	template <typename NormBijection>
	real mutualInformation(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	//! Computes mutual information between two signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation().
	*/

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

}

#include "tim/core/mutual_information.hpp"

#endif
