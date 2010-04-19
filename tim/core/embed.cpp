#include "tim/core/embed.h"
#include "tim/core/signal_tools.h"

namespace Tim
{

	TIM SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);
		
		const integer n = signal->dimension();
		const integer samples = signal->samples();

		const integer embedDimension = k * n;
		const integer embedLag = (k - 1) * dt;
		const integer embedSamples = std::max(samples - embedLag, 0);

		const SignalPtr embedSignal = SignalPtr(
			new Signal(embedSamples, embedDimension, signal->t() + embedLag));

		integer sBegin = 0;
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
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		const integer dimension =
			signal->dimension();
		const integer samples = 
			signal->samples();

		const integer futureShift = 
			delayEmbedFutureShift(k, dt);

		const integer embedSamples = 
			std::max(samples - futureShift, 0);

		const SignalPtr embedSignal = SignalPtr(
			new Signal(embedSamples, dimension,
			signal->t(), &*signal->data().rowBegin(futureShift)));

		/*
		if (embedSamples > 0)
		{
			const integer beginIndex = futureShift * dimension;
			const integer endIndex = beginIndex + embedSamples * dimension;

			std::copy(
				signal->data().rawBegin() + beginIndex,
				signal->data().rawBegin() + endIndex,
				embedSignal->data().rawBegin());
		}
		*/

		return embedSignal;
	}

	TIM integer delayEmbedFutureShift(
		integer k, 
		integer dt)
	{
		PENSURE_OP(k, >, 0);
		PENSURE_OP(dt, >=, 1);

		return dt;
	}

}
