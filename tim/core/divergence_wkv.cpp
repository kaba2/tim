#include "tim/core/divergence_wkv.h"

namespace Tim
{

	TIM real divergenceWkv(
		const Signal& xSignal,
		const Signal& ySignal)
	{
		const Signal xSignalSet[] = {xSignal};
		const Signal ySignalSet[] = {ySignal};

		return Tim::divergenceWkv(
			range(xSignalSet),
			range(ySignalSet));
	}

}
