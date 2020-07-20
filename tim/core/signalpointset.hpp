#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>

namespace Tim
{

	template <ranges::forward_range Signal_Range>
	SignalPointSet::SignalPointSet(
		const Signal_Range& signalSet)
		: kdTree_(Pointer_Locator<dreal>(ranges::empty(signalSet) ? 0 : std::begin(signalSet)->dimension()))
		, pointSet_()
		, signals_(ranges::size(signalSet))
		, samples_(0)
		, windowBegin_(0)
		, windowEnd_(0)
		, dimensionBegin_(0)
		, dimension_(kdTree_.n())
		, timeBegin_(0)
	{
		ENSURE(!ranges::empty(signalSet));
		PENSURE(equalDimension(signalSet));

		createPointSet(signalSet);
	}

	template <ranges::forward_range Signal_Range>
	SignalPointSet::SignalPointSet(
		const Signal_Range& signalSet,
		integer dimensionBegin,
		integer dimensionEnd)
		: kdTree_(Pointer_Locator<dreal>(dimensionEnd - dimensionBegin))
		, pointSet_()
		, signals_(ranges::size(signalSet))
		, samples_(0)
		, windowBegin_(0)
		, windowEnd_(0)
		, dimensionBegin_(dimensionBegin)
		, dimension_(dimensionEnd - dimensionBegin)
		, timeBegin_(0)
	{
		ENSURE(!ranges::empty(signalSet));
		PENSURE(equalDimension(signalSet));
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, std::begin(signalSet)->dimension());

		createPointSet(signalSet);
	}

	// Private

	template <ranges::forward_range Signal_Range>
	void SignalPointSet::createPointSet(
		const Signal_Range& signalSet)
	{
		// Find out the time interval on which
		// all trials are defined.

		Integer2 sharedTime = sharedTimeInterval(signalSet);
		integer tBegin = sharedTime[0];
		integer tEnd = sharedTime[1];
		integer samples = tEnd - tBegin;

		// Store points in an interleaved
		// manner.

		integer signals = ranges::size(signalSet);

		pointSet_.resize(samples * signals);

		auto iter = signalSet.begin();
		for (integer i = 0;i < signals;++i)
		{
			const Signal& signal = *iter;
			for (integer t = tBegin;t < tEnd;++t)
			{
				const dreal* point = 
					std::begin(signal.pointRange(dimensionBegin_))[t - signal.t()];
				pointSet_[(t - signal.t()) * signals + i] = kdTree_.insert(point);
			}
			
			++iter;
		}

		signals_ = signals;
		samples_ = samples;
		timeBegin_ = tBegin;

		windowBegin_ = tBegin;
		windowEnd_ = tBegin + samples;

		kdTree_.refine(SplitRule());
	}

}

#endif
