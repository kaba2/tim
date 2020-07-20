// Description: Merging signals into higher-dimensional signals

#ifndef TIM_SIGNAL_MERGE_H
#define TIM_SIGNAL_MERGE_H

#include "tim/core/signal.h"
#include "tim/core/signal_properties.h"

#include <pastel/sys/range.h>
#include <pastel/sys/array/array.h>

namespace Tim
{

	//! Merges a signal set into a higher-dimensional signal.
	template <
		ranges::forward_range Signal_Range,
		ranges::forward_range Lag_Range>
	SignalData merge(
		const Signal_Range& signalSet,
		const Lag_Range& lagSet)
	{
		ENSURE_OP(ranges::size(signalSet), ==, ranges::size(lagSet));

		if (ranges::empty(signalSet) ||
			ranges::empty(lagSet))
		{
			return SignalData();
		}

		// Compute joint dimension.

		integer jointDimension = 0;
		for (auto&& signal : signalSet) {
			jointDimension += signal.dimension();
		}

		Integer2 sharedTime = 
			sharedTimeInterval(signalSet, lagSet);
		integer samples = sharedTime[1] - sharedTime[0];

		// Allocate the joint signal.

		SignalData jointSignal(samples, jointDimension, sharedTime[0]);
		
		if  (samples == 0)
		{
			// There is no common time interval that
			// all signals would share.
			return jointSignal;
		}

		// Copy the signals into parts of the joint signal.

		integer dimensionOffset = 0;

		auto signalIter = std::begin(signalSet);
		for (integer lag : lagSet) 
		{
			const Signal& signal = *signalIter;
			const integer lagOffset = sharedTime[0] - (signal.t() + lag);
			integer dimension = signal.dimension();

			auto jointSliced = jointSignal.data().slicex(dimensionOffset);

			for (integer i = 0;i < samples;++i)
			{
				ranges::copy(
					signal.data().rowRange(i + lagOffset),
					std::begin(jointSliced.rowRange(i)));
			}

			dimensionOffset += dimension;

			++signalIter;
		}

		return jointSignal;
	}

	//! Merges a signal set into a higher-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(signalSet, constantRange(0, ranges::size(signalSet)));
	See the documentation for that function.
	*/
	template <ranges::forward_range Signal_Range>
	Signal merge(
		const Signal_Range& signalSet)
	{
		return Tim::merge(signalSet,
			constantRange(0, ranges::size(signalSet)));
	}

	//! Merges two signal sets pairwise into a new signal set.
	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Signal_OutputIterator>
	void merge(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		Signal_OutputIterator result,
		integer xLag = 0, integer yLag = 0)
	{
		ENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(ySignalSet));
		
		auto xIter = std::begin(xSignalSet);
		auto xIterEnd = std::end(xSignalSet);
		auto yIter = std::begin(ySignalSet);

		while(xIter != xIterEnd)
		{
			*result = merge(*xIter, *yIter, xLag, yLag);
			
			++result;
			++xIter;
			++yIter;
		}
	}

	//! Merges signals into a single high-dimensional signal.
	/*!
	Preconditions:
	ranges::size(lagSet) == ensembleSet.height()

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
		typename Lag_Range>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result,
		const Lag_Range& lagSet)
	{
		ENSURE_OP(ranges::size(lagSet), ==, ensembleSet.height());

		integer trials = ensembleSet.width();
		for (integer i = 0;i < trials;++i)
		{
			*result = merge(ensembleSet.cColumnRange(i), lagSet);
			++result;
		}
	}

	//! Merges signals into a single high-dimensional signal.
	/*!
	This is a convenience function that calls:
	merge(ensembleSet, result,
		constantRange(0, ensembleSet.height()));
	*/
	template <typename Signal_OutputIterator>
	void merge(
		const Array<Signal>& ensembleSet,
		Signal_OutputIterator result)
	{
		Tim::merge(ensembleSet, result,
			constantRange(0, ensembleSet.height()));
	}

	//! Merges two signals into a higher-dimensional signal.
	inline TIM Signal merge(
		const Signal& xSignal,
		const Signal& ySignal,
		integer xLag = 0,
		integer yLag = 0)
	{
		return Tim::merge(
			range({ xSignal , ySignal }),
			range({ xLag, yLag }));
	}

}

#endif
