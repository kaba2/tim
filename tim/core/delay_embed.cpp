#include "tim/core/delay_embed.h"
#include "tim/core/signal_tools.h"

namespace Tim
{

	TIM Signal delayEmbed(
		const Signal& signal,
		integer k,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);
		
		integer n = signal.dimension();
		integer samples = signal.samples();

		integer embedDimension = k * n;
		integer embedLag = (k - 1) * dt;
		integer embedSamples = std::max(samples - embedLag, (integer)0);

		Signal embedSignal(
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
				std::copy(
					signal.data().cRowBegin(s),
					signal.data().cRowEnd(s),
					embedSignal.data().rowBegin(t) + iBegin);

				s += dt;
				iBegin += n;
			}
			++sBegin;
		}

		return embedSignal;
	}

	TIM Signal delayEmbedFuture(
		const Signal& signal,
		integer dt)
	{
		ENSURE_OP(dt, >=, 1);

		integer dimension =
			signal.dimension();
		integer samples = 
			signal.samples();

		integer embedSamples = 
			std::max(samples - dt, (integer)0);

		Signal embedSignal(
			embedSamples, 
			dimension,
			signal.t(), 
			&*removeConst(signal.data()).rowBegin(dt));

		/*
		if (embedSamples > 0)
		{
			const integer beginIndex = futureShift * dimension;
			const integer endIndex = beginIndex + embedSamples * dimension;

			std::copy(
				signal.data().rawBegin() + beginIndex,
				signal.data().rawBegin() + endIndex,
				embedSignal.data().rawBegin());
		}
		*/

		return embedSignal;
	}

}
