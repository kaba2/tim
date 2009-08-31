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
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			forwardRange(constantIterator(zSignal)),
			yLag,
			zLag,
			kNearest);
	}

}
