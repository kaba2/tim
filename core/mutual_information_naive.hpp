#ifndef TIM_MUTUAL_INFORMATION_NAIVE_HPP
#define TIM_MUTUAL_INFORMATION_NAIVE_HPP

#include "tim/core/differential_entropy.h"

namespace Tim
{

	template <typename NormBijection>
	real mutualInformationFromEntropy(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		std::vector<SignalPtr> marginalSignalSet;
		slice(jointSignal, partition, marginalSignalSet);

		const integer signals = marginalSignalSet.size();
		const integer samples = jointSignal->samples();

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE2(jointSignal->samples() == marginalSignalSet[i]->samples(),
				jointSignal->samples(), marginalSignalSet[i]->samples());

			jointDimension += marginalSignalSet[i]->dimension();
		}

		ENSURE2(jointDimension == jointSignal->dimension(),
			jointDimension, jointSignal->dimension());

		real estimate = 0;

		for (integer i = 0;i < signals;++i)
		{
			estimate += differentialEntropy(
				marginalSignalSet[i],
				kNearest,
				maxRelativeError,
				normBijection);
		}

		estimate -= differentialEntropy(
			jointSignal,
			kNearest,
			maxRelativeError,
			normBijection);

		return estimate;
	}

}

#endif
