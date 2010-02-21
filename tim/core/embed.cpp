#include "tim/core/embed.h"
#include "tim/core/signal_tools.h"

namespace Tim
{

	TIM SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer t0,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(t0, >=, 0);
		ENSURE_OP(dt, >=, 1);

		const integer n = signal->dimension();

		const integer embedDimension = k * n;
		const integer embedSampleWidth = (k - 1) * dt + 1;
		const integer embedSamples = (signal->samples() - t0) - embedSampleWidth + 1;

		const SignalPtr embedSignal = SignalPtr(
			new Signal(embedSamples, embedDimension));

		integer sBegin = t0;
		for (integer t = 0;t < embedSamples;++t)
		{
			integer iBegin = 0;
			integer s = sBegin;
			for (integer j = 0;j < k;++j)
			{
				std::copy(
					signal->data().rowBegin(s),
					signal->data().rowEnd(s),
					embedSignal->data().rowBegin(t) + iBegin);

				s += dt;
				iBegin += n;
			}
			++sBegin;
		}

		return embedSignal;
	}

	TIM SignalPtr delayEmbedFuture(
		const SignalPtr& signal,
		integer k,
		integer t0,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(t0, >=, 0);
		ENSURE_OP(dt, >=, 1);

		const integer futureShift = 
			delayEmbedFutureShift(k, t0, dt);

		return delayEmbed(signal, 1, futureShift);
	}

	TIM integer delayEmbedFutureShift(
		integer k, 
		integer t0, 
		integer dt)
	{
		PENSURE_OP(k, >, 0);
		PENSURE_OP(t0, >=, 0);
		PENSURE_OP(dt, >=, 1);

		return t0 + dt * k;
	}

}
