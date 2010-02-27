#include "tim/core/signal_merge.h"

namespace Tim
{

	TIM SignalPtr merge(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer xLag,
		integer yLag)
	{
		const SignalPtr signalSet[2] = {xSignal, ySignal};
		const integer lagSet[2] = {xLag, yLag};

		return Tim::merge(
			forwardRange(signalSet),
			forwardRange(lagSet));
	}

}
