// Description: Algorithms for Signal's

#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

#include <pastel/gfx/color.h>

#include <pastel/sys/smallset.h>
#include <pastel/sys/forwardrange.h>

#include <iostream>

namespace Tim
{

	TIM std::ostream& operator<<(std::ostream& stream, const Signal& signal);

	//! Returns the minimum number of samples among the signals.
	/*!
	Preconditions:
	SignalPtr_Iterator must dereference to SignalPtr.
	*/
	template <typename SignalPtr_Iterator>
	integer minSamples(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same dimension.
	/*!
	Preconditions:
	SignalPtr_Iterator must dereference to SignalPtr.
	*/
	template <typename SignalPtr_Iterator>
	bool equalDimension(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	// Merge
	// -----

	//! Merges a signal set into a higher-dimensional signal.
	/*!
	Preconditions:
	SignalPtr_Iterator dereferences to SignalPtr.

	The output signal will have as many samples as
	the input signal with the least number of samples.
	*/
	template <
		typename SignalPtr_Iterator,
		typename Integer_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet);

	template <typename SignalPtr_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Merges two signal sets pairwise into a new signal set.
	/*!
	Preconditions:
	SignalPtr_X_Iterator dereferences to SignalPtr.
	SignalPtr_X_Iterator dereferences to SignalPtr.
	SignalPtr_OutputIterator dereferences to SignalPtr.
	yLag >= 0
	*/
	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_OutputIterator>
	void merge(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		SignalPtr_OutputIterator result,
		integer yLag = 0);

	//! Merges signals into a single high-dimensional signal.
	/*!
	Preconditions:
	lagSet.size() == ensembleSet.height()

	ensembleSet:
	An array where each row contains trials of one signal.
	The signals in the same column must correspond in time.

	result:
	An iterator where the merged signal trials are written.

	lagSet:
	For each row of 'ensembleSet', a lag giving the
	delay in samples to apply to the signal in that
	row before the merging.
	*/
	template <
		typename SignalPtr_OutputIterator,
		typename Integer_Iterator>
	void merge(
		const Array<SignalPtr, 2>& ensembleSet,
		SignalPtr_OutputIterator result,
		const ForwardRange<Integer_Iterator> lagSet);

	//! Merges signals into a single high-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(ensembleSet, result,
		forwardRange(constantIterator(0), ensembleSet.height()));
	*/
	template <typename SignalPtr_OutputIterator>
	void merge(
		const Array<SignalPtr, 2>& ensembleSet,
		SignalPtr_OutputIterator result);

	//! Merges two signals into a higher-dimensional signal.
	/*!
	Preconditions:
	SignalPtr_Iterator dereferences to SignalPtr.

	The output signal will have as many samples as
	the input signal with the least number of samples.
	*/
	TIM SignalPtr merge(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer xLag = 0,
		integer yLag = 0);

	// Split
	// -----

	//! Creates aliases for 1d-marginal signals.
	/*!
	See the more specialized
	'split' function for more information.
	*/
	template <typename SignalPtr_OutputIterator>
	void split(
		const SignalPtr& jointSignal,
		SignalPtr_OutputIterator signalSet);

	//! Creates aliases for marginal signals.
	/*!
	The partition of the joint signal is given
	by the 'partition' set. Assume that 'jointSignal'
	is of dimension 4 and 'partition' contains the 
	numbers 0, 2, 3, 4. Then the marginal signals
	are splitd with the dimension subranges
	[0, 2[, [2, 3[, [3, 4[.
	*/
	template <typename SignalPtr_OutputIterator>
	void split(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		SignalPtr_OutputIterator signalSet);

	//! Creates an alias for a marginal signal.
	/*
	Returns a signal S which refers to a marginal signal of P. 
	For each sample x in P, S refers only to the dimension subrange 
	[dimensionBegin, dimensionEnd[. S then has dimension
	(dimensionEnd - dimensionBegin). Note S shares memory
	with P and thus changes in either are reflected in
	the other.
	*/
	TIM SignalPtr split(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd);

	//! Computes the covariance of the signal samples.
	TIM void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result);

	//! Transforms the given signal to identity covariance.
	TIM void normalizeCovariance(
		const SignalPtr& signal,
		const MatrixD& covariance);

	template <typename Image_View>
	void drawSignal(
		const SignalPtr& signal,
		const View<2, Color, Image_View>& image);

}

#include "tim/core/signal_generate.h"

#include "tim/core/signal_tools.hpp"

#endif
