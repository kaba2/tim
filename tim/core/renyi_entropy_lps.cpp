#include "tim/core/renyi_entropy_lps.h"

namespace Tim
{

	TIM real renyiEntropyLps(
		const SignalPtr& signal,
		real q,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::renyiEntropyLps(
			signal, q,
			maxRelativeError,
			kNearest);
	}

}
