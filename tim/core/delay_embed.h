// Description: Delay embedding

#ifndef TIM_DELAY_EMBED_H
#define TIM_DELAY_EMBED_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>

#include <vector>

namespace Tim
{

	//! Delay-embeds a signal into a higher dimensional space.
	/*!
	Preconditions:
	k > 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.
	*/

	inline TIM SignalData delayEmbed(
		const Signal& signal,
		integer k,
		integer dt = 1)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		integer n = signal.dimension();
		integer samples = signal.samples();

		integer embedDimension = k * n;
		integer embedLag = (k - 1) * dt;
		integer embedSamples = std::max(samples - embedLag, (integer)0);

		SignalData embedSignal(
			embedSamples, 
			embedDimension, 
			signal.t() + embedLag);

		integer sBegin = 0;
		for (integer t = 0;t < embedSamples;++t)
		{
			integer iBegin = 0;
			integer s = sBegin;
			for (integer j = 0;j < k;++j)
			{
				ranges::copy(
					signal.data().rowRange(s),
					embedSignal.data().rowRange(t) | ranges::drop(iBegin));

				s += dt;
				iBegin += n;
			}
			++sBegin;
		}

		return embedSignal;
	}


	//! Performs delay embedding for a set of signals.
	/*!
	Preconditions:
	k > 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.
	*/
	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbed(
		Signal_Input inputSet,
		Signal_Output output,
		integer k,
		integer dt = 1)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		while(!inputSet.empty())
		{
			output(delayEmbed(inputSet.get(), k, dt));
			inputSet.pop();
		}
	}

	//! Returns the future of a signal under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	inline TIM Signal delayEmbedFuture(
		const Signal& signal,
		integer dt = 1)
	{
		return Signal(signal.data().slicex(dt), signal.t());
	}

	//! Computes the futures of signals under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbedFuture(
		Signal_Input inputSet,
		Signal_Output output,
		integer dt = 1)
	{
		ENSURE_OP(dt, >=, 1);

		while(!inputSet.empty())
		{
			output(Tim::delayEmbedFuture(inputSet.get(), dt));
			inputSet.pop();
		}
	}

}

#endif
