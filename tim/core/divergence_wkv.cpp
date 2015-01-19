#include "tim/core/divergence_wkv.h"

#include <pastel/sys/input/single_input.h>

namespace Tim
{

	TIM real divergenceWkv(
		const Signal& xSignal,
		const Signal& ySignal)
	{
		return Tim::divergenceWkv(
			range({ &xSignal }),
			range({ &ySignal }));
	}

}
