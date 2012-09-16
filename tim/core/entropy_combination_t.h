// Description: Temporal estimation of entropy combinations

#ifndef TIM_ENTROPY_COMBINATION_T_H
#define TIM_ENTROPY_COMBINATION_T_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/sys/array.h>

namespace Tim
{

	//! Computes a temporal entropy combination of signals.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	odd(filter.size())

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
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	result:
	Temporal estimates of the entropy combination of the
	signals.

	lagSet:
	Lags to apply to each signal.

	filter:
	An array of coefficients by which to weight the results
	in the time-window. The center of the array corresponds 
	to the current time instant. The width of the array can 
	be arbitrary but must be odd. The coefficients must sum
	to a non-zero value.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate entropy combination.

	Returns:
	The temporal estimates in a 1d-signal.
	*/

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator,
		typename Real_Filter_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const boost::iterator_range<Integer_Iterator>& lagSet,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter);

	//! Computes a temporal entropy combination of signals.
	/*!
	This is a convenience function that calls:

	temporalEntropyCombination(
		signalSet,
		rangeSet,
		timeWindowRadius,
		lagSet,
		kNearest,
		constantRange((real)1, 1));

	See the documentation for that function.
	*/

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const boost::iterator_range<Integer_Iterator>& lagSet,
		integer kNearest = 1);

	//! Computes a temporal entropy combination of signals.
	/*!
	This is a convenience function that calls:

	temporalEntropyCombination(
		signalSet,
		rangeSet,
		timeWindowRadius,
		constantRange(0, signalSet.height()));

	See the documentation for that function.
	*/

	template <typename Integer3_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius);

}

#include "tim/core/entropy_combination_t.hpp"

#endif
