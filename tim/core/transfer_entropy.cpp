#include "tim/core/transfer_entropy.h"

namespace Tim
{

	TIM real transferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& wSignal,
		integer xLag,
		integer yLag,
		integer wLag,
		integer kNearest)
	{
		return Tim::transferEntropy(
			constantRange(xSignal),
			constantRange(ySignal),
			constantRange(wSignal),
			xLag, yLag, wLag,
			kNearest);
	}

}
