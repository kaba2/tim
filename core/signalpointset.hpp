#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"

#include <pastel/sys/ensure.h>

namespace Tim
{

	template <typename Signal_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<Signal_Iterator>& signalSet)
		: kdTree_(ofDimension(signalSet.empty() ? 0 : signalSet.front()->dimension()))
		, signalSet_(signalSet.begin(), signalSet.end())
		, pointSet_()
		, objectSet_()
		, samples_(minSamples(signalSet))
		, timeBegin_(samples_)
		, timeEnd_(samples_)
		, dimensionBegin_(0)
		, dimension_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		construct(samples_, samples_, 
			0, signalSet.front()->dimension());
	}

	template <typename Signal_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeBegin,
		integer timeEnd)
		: kdTree_(ofDimension(signalSet.empty() ? 0 : signalSet.front()->dimension()))
		, signalSet_(signalSet.begin(), signalSet.end())
		, pointSet_()
		, objectSet_()
		, samples_(minSamples(signalSet))
		, timeBegin_(samples_)
		, timeEnd_(samples_)
		, dimensionBegin_(0)
		, dimension_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		construct(timeBegin, timeEnd, 
			0, signalSet.front()->dimension());
	}

	template <typename Signal_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeBegin,
		integer timeEnd,
		integer dimensionBegin,
		integer dimensionEnd)
		: kdTree_(ofDimension(dimensionEnd - dimensionBegin))
		, signalSet_(signalSet.begin(), signalSet.end())
		, pointSet_()
		, objectSet_()
		, samples_(minSamples(signalSet))
		, timeBegin_(samples_)
		, timeEnd_(samples_)
		, dimensionBegin_(0)
		, dimension_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		construct(timeBegin, timeEnd,
			dimensionBegin, dimensionEnd);
	}

}

#endif
