#ifndef TIM_SIGNAL_SPLIT_HPP
#define TIM_SIGNAL_SPLIT_HPP

#include "tim/core/signal_split.h"

namespace Tim
{

	template <typename SignalPtr_OutputIterator>
	void split(
		const SignalPtr& jointSignal,
		SignalPtr_OutputIterator signalSet)
	{
		const integer dimension = jointSignal->dimension();

		SmallSet<integer> partition;
		partition.reserve(dimension + 1);
		for (integer i = 0;i <= dimension;++i)
		{
			partition.insert(i);
		}

		Tim::split(jointSignal, partition, signalSet);
	}

	template <typename SignalPtr_OutputIterator>
	void split(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		SignalPtr_OutputIterator signalSet)
	{
		ENSURE_OP(partition.size(), >=, 2);

		const integer dimension = jointSignal->dimension();
		const integer samples = jointSignal->samples();
		const integer signals = partition.size() - 1;

		for (integer x = 0;x < signals;++x)
		{
			const integer marginalDimension = 
				partition[x + 1] - partition[x];

			*signalSet = Tim::split(jointSignal, partition[x], 
				partition[x] + marginalDimension);
			++signalSet;
		}
	}

}

#endif
