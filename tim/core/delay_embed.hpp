#ifndef TIM_DELAY_EMBED_HPP
#define TIM_DELAY_EMBED_HPP

#include "tim/core/delay_embed.h"

namespace Tim
{

	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbed(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
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
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer dt)
	{
		ENSURE_OP(dt, >=, 1);

		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalEnd = signalSet.end();
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = Tim::delayEmbedFuture(*signalIter, dt);
			++outputIter;
			++signalIter;
		}
	}

}

#endif
