#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"

#include <pastel/sys/ensure.h>

namespace Tim
{

	template <typename Signal_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<Signal_Iterator>& signalSet,
		SignalPointSet_TimeWindow::Enum timeWindowStart)
		: kdTree_(ofDimension(signalSet.empty() ? 0 : signalSet.front()->dimension()))
		, signalSet_(signalSet.begin(), signalSet.end())
		, objectSet_()
		, samples_(minSamples(signalSet))
		, timeBegin_(samples_)
		, timeEnd_(samples_)
		, dimensionBegin_(0)
		, dimension_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		construct(timeWindowStart,
			0, signalSet.front()->dimension());
	}

	template <typename Signal_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<Signal_Iterator>& signalSet,
		SignalPointSet_TimeWindow::Enum timeWindowStart,
		integer dimensionBegin,
		integer dimensionEnd)
		: kdTree_(ofDimension(dimensionEnd - dimensionBegin))
		, signalSet_(signalSet.begin(), signalSet.end())
		, objectSet_()
		, samples_(minSamples(signalSet))
		, timeBegin_(samples_)
		, timeEnd_(samples_)
		, dimensionBegin_(0)
		, dimension_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		construct(timeWindowStart,
			dimensionBegin, dimensionEnd);
	}

}

#endif
