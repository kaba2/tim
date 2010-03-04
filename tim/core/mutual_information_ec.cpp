#include "tim/core/mutual_information_ec.h"

namespace Tim
{

	TIM real mutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer xLag,
		integer yLag,
		integer kNearest)
	{
		return Tim::mutualInformation(
			constantRange(xSignal),
			constantRange(ySignal),
			xLag,
			yLag,
			kNearest);
	}

}

