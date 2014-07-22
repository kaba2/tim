#include "tim/core/signal_merge.h"

namespace Tim
{

	TIM Signal merge(
		const Signal& xSignal,
		const Signal& ySignal,
		integer xLag,
		integer yLag)
	{
		const Signal signalSet[2] = {xSignal, ySignal};
		const integer lagSet[2] = {xLag, yLag};

		return Tim::merge(
			range(signalSet),
			range(lagSet));
	}

}
