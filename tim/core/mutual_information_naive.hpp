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
		const boost::iterator_range<Signal_Iterator>& signalSet,
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

		Signal_Iterator iter = signalSet.begin();
		Signal_Iterator iterEnd = signalSet.end();

		real estimate = 0;
		while(iter != iterEnd)
		{
			const Signal signal = *iter;

			estimate += differentialEntropyKl(
				signal,
				kNearest,
				normBijection);

			++iter;
		}

		const Signal jointSignal = merge(signalSet);
		estimate -= differentialEntropyKl(
			jointSignal,
			kNearest,
			normBijection);

		return estimate;
	}

}

#endif
