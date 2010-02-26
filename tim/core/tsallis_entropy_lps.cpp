#include "tim/core/tsallis_entropy_lps.h"

namespace Tim
{

	TIM real tsallisEntropyLps(
		const SignalPtr& signal,
		real q,
		integer kNearest)
	{
		return Tim::tsallisEntropyLps(
			signal, q,
			kNearest);
	}

	TIM real tsallisDecideK(real q, integer kNearestSuggestion)
	{
		PENSURE_OP(q, >, 0);
		PENSURE_OP(kNearestSuggestion, >=, 0);

		integer kNearest = kNearestSuggestion;

		if (kNearestSuggestion == 0)
		{
			// We get to decide the k.

			kNearest = 2 * std::ceil(q);
		}
		else if (kNearestSuggestion <= q - 1)
		{
			// The algorithm is not defined
			// for such k. Find the smallest
			// k for which the algorithm is
			// defined.

			if (kNearestSuggestion < q - 1)
			{
				// 0 < kNearestSuggestion
				// thus
				// 0 < q - 1

				kNearest = std::ceil(q - 1);
			}
			else
			{
				kNearest = kNearestSuggestion + 1;
			}
		}

		return kNearest;
	}

}
