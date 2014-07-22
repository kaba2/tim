#ifndef TIM_MUTUAL_INFORMATION_NAIVE_HPP
#define TIM_MUTUAL_INFORMATION_NAIVE_HPP

#include "tim/core/differential_entropy.h"

namespace Tim
{

	template <
		typename SignalPtr_Range,
		typename NormBijection,
		typename Real_OutputIterator>
	real mutualInformationFromEntropy(
		const SignalPtr_Range& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection,
		Real_OutputIterator result)
	{
		ENSURE_OP(kNearest, >, 0);

		if (signalSet.empty())
		{
			return 0;
		}

		auto iter = signalSet.begin();
		auto iterEnd = signalSet.end();

		real estimate = 0;
		while(iter != iterEnd)
		{
			const Signal& signal = *iter;

			estimate += differentialEntropyKl(
				signal,
				kNearest,
				normBijection);

			++iter;
		}

		const Signal& jointSignal = merge(signalSet);
		estimate -= differentialEntropyKl(
			jointSignal,
			kNearest,
			normBijection);

		return estimate;
	}

}

#endif
