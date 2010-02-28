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
		
		// We want the length of the delay-embedded 
		// signal to be independent of k. To do this
		// we extend the finite-length signal to an
		// infinite length one by repetition, i.e. turn 
		// the signal into a periodic one.

		const integer n = signal->dimension();
		const integer samples = signal->samples();

		const integer embedDimension = k * n;
		const integer embedSamples = samples - t0;
		//const integer embedSampleWidth = (k - 1) * dt + 1;
		//const integer embedSamples = (signal->samples() - t0) - embedSampleWidth + 1;

		if (embedSamples <= 0)
		{
			// The embedding shift goes out of the
			// signal. Return an empty signal.
			return SignalPtr(new Signal(0, embedDimension));
		}

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
				if (s > samples)
				{
					s %= samples;
				}
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
