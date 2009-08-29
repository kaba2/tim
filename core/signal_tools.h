// Description: Tools for working with Signal's

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

	TIMCORE std::ostream& operator<<(std::ostream& stream, const Signal& signal);

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

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_OutputIterator>
	void merge(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		SignalPtr_OutputIterator result,
		integer yLag = 0,
		integer zLag = 0);

	template <typename SignalPtr_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Merges two signals into a higher-dimensional signal.
	/*!
	Preconditions:
	SignalPtr_Iterator dereferences to SignalPtr.

	The output signal will have as many samples as
	the input signal with the least number of samples.
	*/
	TIMCORE SignalPtr merge(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer yLag = 0);

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
	TIMCORE SignalPtr split(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd);

	//! Computes the covariance of the signal samples.
	TIMCORE void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result);

	//! Transforms the given signal to identity covariance.
	TIMCORE void normalizeCovariance(
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
