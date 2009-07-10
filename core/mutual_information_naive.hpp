#ifndef TIM_MUTUAL_INFORMATION_NAIVE_HPP
#define TIM_MUTUAL_INFORMATION_NAIVE_HPP

#include "tim/core/differential_entropy.h"

namespace Tim
{

	template <typename NormBijection>
	real mutualInformationFromEntropy(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		const integer signals = signalSet.size();
		const SignalPtr jointSignal = merge(signalSet);

		real estimate = 0;

		for (integer i = 0;i < signals;++i)
		{
			estimate += differentialEntropy(
				signalSet[i],
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
