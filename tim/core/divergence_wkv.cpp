#include "tim/core/divergence_wkv.h"

namespace Tim
{

	TIM real divergenceWkv(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal)
	{
		const SignalPtr xSignalSet[] = {xSignal};
		const SignalPtr ySignalSet[] = {ySignal};

		return Tim::divergenceWkv(
			forwardRange(xSignalSet),
			forwardRange(ySignalSet));
	}

}
