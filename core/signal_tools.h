#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

#include <pastel/sys/smallset.h>

#include <iostream>

namespace Tim
{

	TIMCORE std::ostream& operator<<(std::ostream& stream, const Signal& signal);

	//! Merges signals into a higher-dimensional signal.
	/*!
	Note that in contrast to slicing, merging needs 
	to allocate new memory and copy the signals.
	The output signal will have as many samples as
	the input signal with the least number of samples.
	*/
	TIMCORE SignalPtr merge(
		const std::vector<SignalPtr>& signalList);

	TIMCORE SignalPtr merge(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal);

	//! Creates aliases for 1d-marginal signals.
	/*!
	See the more specialized
	'slice' function for more information.
	*/
	TIMCORE void slice(
		const SignalPtr& jointSignal,
		std::vector<SignalPtr>& marginalSet);

	//! Creates aliases for marginal signals.
	/*!
	The partition of the joint signal is given
	by the 'partition' set. Assume that 'jointSignal'
	is of dimension 4 and 'partition' contains the 
	numbers 0, 2, 3, 4. Then the marginal signals
	are sliced with the dimension subranges
	[0, 2[, [2, 3[, [3, 4[, see the more specialized
	'slice' function for more information.
	*/
	TIMCORE void slice(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		std::vector<SignalPtr>& marginalSet);

	//! Creates an alias for a marginal signal.
	/*
	Returns a signal S which refers to a marginal signal of P. 
	For each sample x in P, S refers only to the dimension subrange 
	[dimensionBegin, dimensionEnd[. S then has dimension
	(dimensionEnd - dimensionBegin). Note S shares memory
	with P and thus changes in either are reflected in
	the other.
	*/
	TIMCORE SignalPtr slice(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd);

	//! Creates a point set of the signal samples.
	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& pointSet);

	//! Creates a point set of the signal samples.
	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<PointD>& pointSet);

	//! Creates a point set of the signal samples.
	TIMCORE void constructPointSet(
		const std::vector<SignalPtr>& ensemble,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<PointD>& pointSet);

	//! Computes the covariance of the signal samples.
	TIMCORE void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result);

	//! Transforms the given signal to identity covariance.
	TIMCORE void normalizeCovariance(
		const SignalPtr& signal,
		const MatrixD& covariance);

}

#include "tim/core/signal_generate.h"

#endif
