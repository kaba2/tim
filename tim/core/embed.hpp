#ifndef TIM_EMBED_HPP
#define TIM_EMBED_HPP

#include "tim/core/embed.h"

namespace Tim
{

	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbed(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalEnd = signalSet.end();
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = delayEmbed(*signalIter, k, dt);
			++outputIter;
			++signalIter;
		}
	}

	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbedFuture(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalEnd = signalSet.end();
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = Tim::delayEmbedFuture(*signalIter, k, dt);
			++outputIter;
			++signalIter;
		}
	}

}

#endif
