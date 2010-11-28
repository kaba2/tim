#include "tim/core/signal_split.h"

namespace Tim
{

	TIM SignalPtr split(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signal->dimension());

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = signal->samples();

		SignalPtr splitSignal(new Signal(
			signal->samples(), dimension));

		splitSignal->data() = signal->data()(
			Vector2i(0, dimensionBegin), Vector2i(samples, dimensionEnd));

		return splitSignal;
	}

}