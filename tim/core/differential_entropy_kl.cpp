#include "tim/core/differential_entropy_kl.h"

namespace Tim
{

	TIM real differentialEntropyKl(
		const SignalPtr& signal,
		integer kNearest)
	{
		return Tim::differentialEntropyKl(
			signal, 
			kNearest, Default_NormBijection());
	}

}
