#ifndef TIM_SIGNAL_MERGE_H
#define TIM_SIGNAL_MERGE_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>
#include <pastel/sys/array.h>

namespace Tim
{

	//! Merges a signal set into a higher-dimensional signal.
	/*!
	Preconditions:
	SignalPtr_Iterator dereferences to SignalPtr.

	tMin:
	If given, is filled with the start of the
	time interval from which the merged signal
	was formed from.

	The output signal will have as many samples as
	the input signal with the least number of samples.
	*/
	template <
		typename SignalPtr_Iterator,
		typename Integer_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer* tMin = 0);

	/*!
	This is a convenience function that calls:

	merge(signalSet,
		forwardRange(constantIterator(0)),
		tMin);

	See the documentation for that function.
	*/
	template <typename SignalPtr_Iterator>
	SignalPtr merge(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer* tMin = 0);

	//! Merges two signal sets pairwise into a new signal set.
	/*!
	Preconditions:
	SignalPtr_X_Iterator dereferences to SignalPtr.
	SignalPtr_X_Iterator dereferences to SignalPtr.
	SignalPtr_OutputIterator dereferences to SignalPtr.
	*/
	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_OutputIterator>
	void merge(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		SignalPtr_OutputIterator result,
		integer xLag = 0,
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
		const ForwardRange<Integer_Iterator>& lagSet);

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

}

#include "tim/core/signal_merge.hpp"

#endif
