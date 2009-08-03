#ifndef TIM_EMBED_HPP
#define TIM_EMBED_HPP

#include "tim/core/embed.h"

namespace Tim
{

	template <typename Signal_ForwardIterator, typename OutputIterator>
	void delayEmbed(
		const Signal_ForwardIterator& signalBegin,
		const Signal_ForwardIterator& signalEnd,
		const OutputIterator& outputBegin,
		integer k,
		integer shift,
		integer step)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(shift, >=, 0);
		ENSURE_OP(step, >=, 1);

		Signal_ForwardIterator signalIter = signalBegin;
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = delayEmbed(*signalIter, k, shift, step);
			++outputIter;
			++signalIter;
		}
	}

	template <typename Signal_ForwardIterator, typename OutputIterator>
	void delayEmbedFuture(
		const Signal_ForwardIterator& signalBegin,
		const Signal_ForwardIterator& signalEnd,
		const OutputIterator& outputBegin,
		integer k,
		integer shift,
		integer step)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(shift, >=, 0);
		ENSURE_OP(step, >=, 1);

		Signal_ForwardIterator signalIter = signalBegin;
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = Tim::delayEmbedFuture(*signalIter, k, shift, step);
			++outputIter;
			++signalIter;
		}
	}

}

#endif
