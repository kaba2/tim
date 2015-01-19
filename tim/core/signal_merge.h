// Description: Merging signals into higher-dimensional signals

#ifndef TIM_SIGNAL_MERGE_H
#define TIM_SIGNAL_MERGE_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/sys/array/array.h>

namespace Tim
{

	//! Merges a signal set into a higher-dimensional signal.
	template <
		typename SignalPtr_Range,
		typename Integer_Iterator>
	Signal merge(
		const SignalPtr_Range& signalSet,
		const boost::iterator_range<Integer_Iterator>& lagSet);

	//! Merges a signal set into a higher-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(signalSet, constantRange(0, signalSet.size()));
	See the documentation for that function.
	*/
	template <typename SignalPtr_Range>
	Signal merge(
		const SignalPtr_Range& signalSet);

	//! Merges two signal sets pairwise into a new signal set.
	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Signal_OutputIterator>
	void merge(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		Signal_OutputIterator result,
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
		typename Signal_OutputIterator,
		typename Integer_Iterator>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result,
		const boost::iterator_range<Integer_Iterator>& lagSet);

	//! Merges signals into a single high-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(ensembleSet, result,
		constantRange(0, ensembleSet.height()));
	*/
	template <typename Signal_OutputIterator>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result);

	//! Merges two signals into a higher-dimensional signal.
	TIM Signal merge(
		const Signal& xSignal,
		const Signal& ySignal,
		integer xLag = 0,
		integer yLag = 0);

}

#include "tim/core/signal_merge.hpp"

#endif
