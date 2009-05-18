#ifndef TIM_MUTUAL_INFORMATION_H
#define TIM_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

namespace Tim
{

	//! Computes mutual information between a set of signals.

	real mutualInformation(
		const std::vector<SignalPtr>& signalSet);

	//! Computes mutual information between two signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation().
	*/

	TIMCORE real mutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal);

}

#include "tim/core/mutual_information.hpp"

#endif
