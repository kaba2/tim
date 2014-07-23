#ifndef TIM_SIGNAL_SPLIT_HPP
#define TIM_SIGNAL_SPLIT_HPP

#include "tim/core/signal_split.h"

namespace Tim
{

	template <typename Signal_OutputIterator>
	void split(
		const Signal& jointSignal,
		Signal_OutputIterator signalSet)
	{
		integer dimension = jointSignal.dimension();

		std::vector<integer> partition;
		partition.reserve(dimension + 1);
		for (integer i = 0;i <= dimension;++i)
		{
			partition.push_back(i);
		}

		Tim::split(jointSignal, partition, signalSet);
	}

	template <typename Signal_OutputIterator>
	void split(

		const Signal& jointSignal,
		const std::vector<integer>& partition,
		Signal_OutputIterator signalSet)
	{
		ENSURE_OP(partition.size(), >=, 2);

		integer signals = partition.size() - 1;

		for (integer x = 0;x < signals;++x)
		{
			PENSURE_OP(partition[x], <, partition[x + 1]);

			integer marginalDimension = 
				partition[x + 1] - partition[x];


			*signalSet = Tim::split(jointSignal, partition[x], 
				partition[x] + marginalDimension);
			++signalSet;
		}
	}

}

#endif
