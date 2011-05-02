#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/ensure.h>

namespace Tim
{

	template <typename SignalPtr_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet)
		: kdTree_(Array_PointPolicy<real>(signalSet.empty() ? 0 : signalSet.front()->dimension()))
		, pointSet_()
		, signals_(signalSet.size())
		, samples_(0)
		, windowBegin_(0)
		, windowEnd_(0)
		, dimensionBegin_(0)
		, dimension_(kdTree_.dimension())
		, timeBegin_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		createPointSet(signalSet);
	}

	template <typename SignalPtr_Iterator>
	SignalPointSet::SignalPointSet(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer dimensionBegin,
		integer dimensionEnd)
		: kdTree_(Array_PointPolicy<real>(dimensionEnd - dimensionBegin))
		, pointSet_()
		, signals_(signalSet.size())
		, samples_(0)
		, windowBegin_(0)
		, windowEnd_(0)
		, dimensionBegin_(dimensionBegin)
		, dimension_(dimensionEnd - dimensionBegin)
		, timeBegin_(0)
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signalSet.front()->dimension());

		createPointSet(signalSet);
	}

	// Private

	template <typename SignalPtr_Iterator>
	void SignalPointSet::createPointSet(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet)
	{
		// Find out the time interval on which
		// all trials are defined.

		const Integer2 sharedTime = sharedTimeInterval(signalSet);
		const integer tBegin = sharedTime[0];
		const integer tEnd = sharedTime[1];
		const integer samples = tEnd - tBegin;

		// Store points in an interleaved
		// manner.

		const integer signals = signalSet.size();
		pointSet_.resize(samples * signals);

		SignalPtr_Iterator iter = signalSet.begin();
		for (integer i = 0;i < signals;++i)
		{
			SignalPtr signal = *iter;
			for (integer t = tBegin;t < tEnd;++t)
			{
				const real* point = 
					signal->pointBegin(dimensionBegin_)[t - signal->t()];
				pointSet_[(t - signal->t()) * signals + i] = kdTree_.insert(point);
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
