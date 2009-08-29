#include "tim/core/differential_entropy_kl.h"

namespace Tim
{

	TIM real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::differentialEntropy(
			signal, maxRelativeError,
			kNearest, Euclidean_NormBijection<real>());
	}

}
