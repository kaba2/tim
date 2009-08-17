#ifndef TIM_MUTUAL_INFORMATION_NAIVE_HPP
#define TIM_MUTUAL_INFORMATION_NAIVE_HPP

#include "tim/core/differential_entropy.h"

namespace Tim
{

	template <
		typename Signal_Iterator,
		typename NormBijection,
		typename Real_OutputIterator>
	real mutualInformationFromEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection,
		Real_OutputIterator result)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		if (signalSet.empty())
		{
			return 0;
		}

		Signal_Iterator iter = signalSet.begin();
		Signal_Iterator iterEnd = signalSet.end();

		std::vector<real> estimate;

		real estimate = 0;
		while(iter != iterEnd)
		{
			const SignalPtr signal = *iter;

			estimate += differentialEntropy(
				signal,
				maxRelativeError,
				kNearest,
				normBijection);

			++iter;
		}

		const SignalPtr jointSignal = merge(signalSet);
		estimate -= differentialEntropy(
			jointSignal,
			maxRelativeError,
			kNearest,
			normBijection);

		return estimate;
	}

}

#endif
