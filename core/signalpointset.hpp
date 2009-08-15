#ifndef TIM_SIGNALPOINTSET_HPP
#define TIM_SIGNALPOINTSET_HPP

#include "tim/core/signalpointset.h"

#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>

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
	{
		ENSURE(!signalSet.empty());
		PENSURE(equalDimension(signalSet));

		// Insert all the points into the tree.

		const integer signals = signalSet.size();
		for (integer t = 0;t < samples_;++t)
		{
			for (integer i = 0;i < signals;++i)
			{
				if (timeWindowStart == SignalPointSet_TimeWindow::StartEmpty)
				{
					kdTree_.insert(
						signalSet_[i]->pointBegin()[t]);
				}
				else
				{
					objectSet_.push_back(
						kdTree_.insert(
						signalSet_[i]->pointBegin()[t]));
				}
			}
		}

		// Compute a fine subdivision for the points.

		kdTree_.refine(SlidingMidpoint2_SplitRule_PointKdTree());

		if (timeWindowStart == SignalPointSet_TimeWindow::StartEmpty)
		{
			// Remove all objects but leave subdivision intact.

			kdTree_.eraseObjects();
		}
		else
		{
			timeBegin_ = 0;
			timeEnd_ = samples_;
		}
	}

}

#endif
