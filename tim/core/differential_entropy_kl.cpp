#include "tim/core/differential_entropy_kl.h"

namespace Tim
{

	TIM real differentialEntropyKl(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::differentialEntropyKl(
			signal, maxRelativeError,
			kNearest, Default_NormBijection());
	}

}
