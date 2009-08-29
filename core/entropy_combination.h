// Description: Estimation of entropy combinations
// Detail: Both temporal and non-temporal variants

#ifndef TIM_ENTROPY_COMBINATION_H
#define TIM_ENTROPY_COMBINATION_H

#include "tim/core/mytypes.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Compute temporal entropy combination.

	template <
		typename Signal_Iterator,
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	void temporalEntropyCombination(
		const ForwardRange<Signal_Iterator>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest = 1);

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	void temporalEntropyCombination(
		const SignalPtr& signal,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest = 1);

	//! Compute entropy combination.

	template <
		typename Signal_Iterator,
		typename Integer3_Iterator>
	real entropyCombination(
		const ForwardRange<Signal_Iterator>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer kNearest = 1);

	template <typename Integer3_Iterator>
	real entropyCombination(
		const SignalPtr& signal,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer kNearest = 1);

}

#include "tim/core/entropy_combination.hpp"

#endif
