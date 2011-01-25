#ifndef TIM_SIGNAL_MERGE_H
#define TIM_SIGNAL_MERGE_H

#include "tim/core/signal.h"

#include <pastel/sys/iterator_range.h>
#include <pastel/sys/array.h>

namespace Tim
{

	//! Merges a signal set into a higher-dimensional signal.
	template <
		typename SignalPtr_Iterator,
		typename Integer_Iterator>
	SignalPtr merge(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		const ForwardIterator_Range<Integer_Iterator>& lagSet);

	//! Merges a signal set into a higher-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(signalSet, constantRange(0, signalSet.size()));
	See the documentation for that function.
	*/
	template <typename SignalPtr_Iterator>
	SignalPtr merge(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet);

	//! Merges two signal sets pairwise into a new signal set.
	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_OutputIterator>
	void merge(
		const ForwardIterator_Range<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardIterator_Range<SignalPtr_Y_Iterator>& ySignalSet,
		SignalPtr_OutputIterator result,
		integer xLag = 0, integer yLag = 0);

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
		const Array<SignalPtr>& ensembleSet,
		SignalPtr_OutputIterator result,
		const ForwardIterator_Range<Integer_Iterator>& lagSet);

	//! Merges signals into a single high-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(ensembleSet, result,
		constantRange(0, ensembleSet.height()));
	*/
	template <typename SignalPtr_OutputIterator>
	void merge(
		const Array<SignalPtr>& ensembleSet,
		SignalPtr_OutputIterator result);

	//! Merges two signals into a higher-dimensional signal.
	TIM SignalPtr merge(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer xLag = 0,
		integer yLag = 0);

}

#include "tim/core/signal_merge.hpp"

#endif
