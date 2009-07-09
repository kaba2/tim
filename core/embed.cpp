#include "tim/core/embed.h"
#include "tim/core/signal_tools.h"

namespace Tim
{

	TIMCORE SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer shift,
		integer step)
	{
		ENSURE1(k > 0, k);
		ENSURE1(shift >= 0, shift);
		ENSURE1(step >= 1, step);

		// Given is a signal S : Z -> R^n.
		// Form a signal R : Z -> R^n : 
		// R(t) = (S(t0 + t), S(t0 + t + dt), ..., S(t0 + t + dt * (k - 1)).
		//
		// Then R is the delay-embedding and
		// t0 is the embedding shift
		// dt is the embedding delay
		// k is the 'embedding factor'
		// d = k n is the embedding dimension

		const integer n = signal->dimension();

		const integer embedDimension = k * n;
		const integer embedSampleWidth = (k - 1) * step + 1;
		const integer embedSamples = (signal->samples() - shift) - embedSampleWidth + 1;

		const SignalPtr embedSignal = SignalPtr(
			new Signal(embedDimension, embedSamples));

		integer sBegin = shift;
		for (integer t = 0;t < embedSamples;++t)
		{
			integer iBegin = 0;
			integer s = sBegin;
			for (integer j = 0;j < k;++j)
			{
				std::copy(
					signal->data().columnBegin(s),
					signal->data().columnEnd(s),
					embedSignal->data().columnBegin(t) + iBegin);

				s += step;
				iBegin += n;
			}
			++sBegin;
		}

		return embedSignal;
	}

	/*
	TIMCORE SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer shift,
		integer step)
	{
		ENSURE1(k > 0, k);
		ENSURE1(shift >= 0, shift);
		ENSURE1(step >= 1, step);

		// Given is a signal S : Z -> R^n.
		// Form a signal R : Z -> R^n : R(t) = S(t0 + dt * t).
		// Form a signal Q : Z -> R^d : Q(t) = (R(t), R(t + 1), ..., R(t + k - 1)).
		//
		// Then Q is the delay-embedding and
		// t0 is the embedding shift
		// dt is the embedding delay
		// k is the 'embedding factor'
		// d = k n is the embedding dimension

		const integer n = signal->dimension();

		const integer embedDimension = k * n;
		const integer culledSamples =
			((signal->samples() - shift) + (step - 1)) / step;
		const integer embedSamples = culledSamples - k + 1;

		const SignalPtr embedSignal = SignalPtr(
			new Signal(embedDimension, embedSamples));

		integer sBegin = shift;
		for (integer t = 0;t < embedSamples;++t)
		{
			integer iBegin = 0;
			integer s = sBegin;
			for (integer j = 0;j < k;++j)
			{
				std::copy(
					signal->data().columnBegin(s),
					signal->data().columnEnd(s),
					embedSignal->data().columnBegin(t) + iBegin);

				s += step;
				iBegin += n;
			}
			sBegin += step;
		}

		return embedSignal;
	}
	*/

}
