// Description: Differential entropy estimation
// Detail: Kozachenko-Leonenko nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_H
#define TIM_DIFFERENTIAL_ENTROPY_KL_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/sys/iterator/constant_iterator.h>

#include <pastel/math/normbijection/normbijection_concept.h>

namespace Tim
{

	//! Temporal differential entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	The norm to use.
	*/
	template <
		typename SignalPtr_Range, 
		typename NormBijection = Default_NormBijection,
		typename Real_Range = ConstantRange<real>>
	Signal temporalDifferentialEntropyKl(
		const SignalPtr_Range& signalSet,
		integer timeWindowRadius,
		integer kNearest = 1,
		const NormBijection& normBijection = NormBijection(),
		const Real_Range& filter = constantRange((real)1, 1));

	//! Differential entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	signalSet contains Signal's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	normBijection:
	The norm to use.

	Returns:
	A differential entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/
	template <
		typename SignalPtr_Range, 
		typename NormBijection = Default_NormBijection>
	real differentialEntropyKl(
		const SignalPtr_Range& signalSet,
		integer kNearest = 1,
		const NormBijection& normBijection = NormBijection());

}

#include "tim/core/differential_entropy_kl.hpp"

#endif
