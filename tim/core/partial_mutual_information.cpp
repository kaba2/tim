#include "tim/core/partial_mutual_information.h"

namespace Tim
{

	TIM real partialMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		return Tim::partialMutualInformation(
			constantRange(xSignal),
			constantRange(ySignal),
			constantRange(zSignal),
			yLag,
			zLag,
			kNearest);
	}

}
