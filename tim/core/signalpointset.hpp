#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>

namespace Tim
{

	template <typename SignalPtr_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		bool startFull)
		: kdTree_(ofDimension(signalSet.empty() ? 0 : signalSet.front()->dimension()))
		, pointSet_()
		, objectSet_()
		, signals_(signalSet.size())
		, samples_(0)
		, timeBegin_(0)
		, timeEnd_(0)
		, dimensionBegin_(0)
		, dimension_(kdTree_.dimension())
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		createPointSet(signalSet);

		construct(startFull);
	}

	template <typename SignalPtr_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		bool startFull,
		integer dimensionBegin,
		integer dimensionEnd)
		: kdTree_(ofDimension(dimensionEnd - dimensionBegin))
		, pointSet_()
		, objectSet_()
		, signals_(signalSet.size())
		, samples_(0)
		, timeBegin_(0)
		, timeEnd_(0)
		, dimensionBegin_(dimensionBegin)
		, dimension_(dimensionEnd - dimensionBegin)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signalSet.front()->dimension());

		createPointSet(signalSet);

		construct(startFull);
	}

	// Private

	template <typename SignalPtr_Iterator>
	void SignalPointSet::createPointSet(
		const ForwardRange<SignalPtr_Iterator>& signalSet)
	{
		const integer signals = signalSet.size();
		const integer samples = minSamples(signalSet);

		// Store points in an interleaved
		// manner.

		pointSet_.resize(samples * signals);

		SignalPtr_Iterator iter = signalSet.begin();
		for (integer i = 0;i < signals;++i)
		{
			SignalPtr signal = *iter;
			for (integer t = 0;t < samples;++t)
			{
				pointSet_[t * signals + i] = 
					signal->pointBegin(dimensionBegin_)[t];
			}
			
			++iter;
		}

		signals_ = signals;
		samples_ = samples;
	}

}

#endif
