#include "tim/core/divergence_wkv.h"

namespace Tim
{

	TIM dreal divergenceWkv(
		const Signal& xSignal,
		const Signal& ySignal)
	{
		return Tim::divergenceWkv(
			range({ xSignal }),
			range({ ySignal }));
	}

}
