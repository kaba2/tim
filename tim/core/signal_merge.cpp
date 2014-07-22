#include "tim/core/signal_merge.h"

namespace Tim
{

	TIM Signal merge(
		const Signal& xSignal,
		const Signal& ySignal,
		integer xLag,
		integer yLag)
	{
		return Tim::merge(
			range({ &xSignal , &ySignal }),
			range({ xLag, yLag }));
	}

}
