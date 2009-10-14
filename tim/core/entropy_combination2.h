// Description: Estimation of entropy combinations
// Detail: Both temporal and non-temporal variants

#ifndef TIM_ENTROPY_COMBINATION2_H
#define TIM_ENTROPY_COMBINATION2_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>
#include <pastel/sys/array.h>

namespace Tim
{

	//! Computes an entropy combination of signals.
	/*!
	Preconditions:
	kNearest > 0

	signalSet:
	An ensemble of joint signals representing trials
	of the same experiment. Note: all the marginal signals
	share the memory with these joint signals.

	rangeSet:
	A sequence of m triples T_i = (a_i, b_i, s_i), 
	where [a_i, b_i] is an interval such that picking those 
	dimensions from the joint signal X gives the marginal 
	signal X_i. The s_i is the factor by which the differential 
	entropy of such a marginal signal is multiplied before summing
	to the end-result.

	timeWindowRadius:
	The radius of the time window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate entropy combination.

	Returns:
	An estimate of the entropy combination of the signals
	C(X_1, ..., X_m) = sum_{i = 1}^m s_i H(X_i) - H(X),
	where H(X) is the differential entropy of X.
	*/

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	real entropyCombination2(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer kNearest = 1);

	//! Computes an entropy combination of signals.
	/*!
	This is a convenience function that calls:

	entropyCombination(
		forwardRange(constantIterator(signal)),
		rangeSet,
		forwardRange(constantIterator(0), signalSet.height()),
		kNearest);

	See the documentation for that function.
	*/

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	real entropyCombination2(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet);

}

#include "tim/core/entropy_combination2.hpp"

#endif
