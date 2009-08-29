#include "tim/core/mutual_information_ec.h"

namespace Tim
{

	TIM real mutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer yLag,
		integer kNearest)
	{
		return Tim::mutualInformation(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			yLag,
			kNearest);
	}

}

