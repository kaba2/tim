#ifndef TIM_MUTUAL_INFORMATION_NAIVE_HPP
#define TIM_MUTUAL_INFORMATION_NAIVE_HPP

#include "tim/core/differential_entropy.h"

namespace Tim
{

	template <
		typename SignalPtr_Iterator,
		typename NormBijection,
		typename Real_OutputIterator>
	real mutualInformationFromEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
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

		SignalPtr_Iterator iter = signalSet.begin();
		SignalPtr_Iterator iterEnd = signalSet.end();

		std::vector<real> estimate;

		real estimate = 0;
		while(iter != iterEnd)
		{
			const SignalPtr signal = *iter;

			estimate += differentialEntropyKl(
				signal,
				kNearest,
				normBijection);

			++iter;
		}

		const SignalPtr jointSignal = merge(signalSet);
		estimate -= differentialEntropyKl(
			jointSignal,
			kNearest,
			normBijection);

		return estimate;
	}

}

#endif
