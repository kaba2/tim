// Description: Algorithms for Signal's

#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"
#include "tim/core/signal_merge.h"
#include "tim/core/signal_properties.h"

#include <pastel/math/matrix.h>

#include <pastel/gfx/color.h>

#include <pastel/sys/smallset.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/array.h>

#include <iostream>

namespace Tim
{

	TIM std::ostream& operator<<(std::ostream& stream, const Signal& signal);

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

	template <typename SignalPtr_Iterator>
	void constructPointSet(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<const real*>& pointSet);

}

#include "tim/core/signal_generate.h"

#include "tim/core/signal_tools.hpp"

#endif
