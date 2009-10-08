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
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			forwardRange(constantIterator(wSignal)),
			xLag, yLag, wLag,
			kNearest);
	}

}
