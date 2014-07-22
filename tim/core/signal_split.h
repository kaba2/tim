// Description: Splitting a signal into lower-dimensional signals

#ifndef TIM_SIGNAL_SPLIT_H
#define TIM_SIGNAL_SPLIT_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

namespace Tim
{

	//! Creates aliases for 1d-marginal signals.
	/*!
	See the more specialized
	'split' function for more information.
	*/
	template <typename Signal_OutputIterator>
	void split(
		const Signal& jointSignal,
		Signal_OutputIterator signalSet);

	//! Creates aliases for marginal signals.
	/*!
	Preconditions:
	partition[x] < partition[x + 1]

	The partition of the joint signal is given
	by the 'partition' set. Assume that 'jointSignal'
	is of dimension 4 and 'partition' contains the 
	numbers 0, 2, 3, 4. Then the marginal signals
	are split with the dimension subranges
	[0, 2[, [2, 3[, [3, 4[.
	*/
	template <typename Signal_OutputIterator>
	void split(
		const Signal& jointSignal,
		const std::vector<integer>& partition,
		Signal_OutputIterator signalSet);

	//! Creates an alias for a marginal signal.
	/*
	Returns a signal S which refers to a marginal signal of P. 
	For each sample x in P, S refers only to the dimension subrange 
	[dimensionBegin, dimensionEnd[. S then has dimension
	(dimensionEnd - dimensionBegin). Note S shares memory
	with P and thus changes in either are reflected in
	the other.
	*/
	TIM Signal split(
		const Signal& signal,
		integer dimensionBegin,
		integer dimensionEnd);

}

#include "tim/core/signal_split.hpp"

#endif
