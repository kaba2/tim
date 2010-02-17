#include "tim/core/tsallis_entropy_lps.h"

namespace Tim
{

	TIM real tsallisEntropyLps(
		const SignalPtr& signal,
		real q,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::tsallisEntropyLps(
			signal, q,
			maxRelativeError,
			kNearest);
	}

}
