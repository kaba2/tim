#include "tim/core/differential_entropy_kl.h"

namespace Tim
{

	TIMCORE real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::differentialEntropy(
			signal, maxRelativeError,
			kNearest, Euclidean_NormBijection<real>());
	}

}
