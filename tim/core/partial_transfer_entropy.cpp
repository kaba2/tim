#include "tim/core/partial_transfer_entropy.h"

namespace Tim
{

	TIM real partialTransferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		const SignalPtr& wSignal,
		integer xLag,
		integer yLag,
		integer zLag,
		integer wLag,
		integer kNearest)
	{
		return Tim::partialTransferEntropy(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			forwardRange(constantIterator(zSignal)),
			forwardRange(constantIterator(wSignal)),
			xLag, yLag, zLag, wLag,
			kNearest);
	}

}
