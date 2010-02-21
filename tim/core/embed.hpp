#ifndef TIM_EMBED_HPP
#define TIM_EMBED_HPP

#include "tim/core/embed.h"

namespace Tim
{

	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbed(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer t0,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(t0, >=, 0);
		ENSURE_OP(dt, >=, 1);

		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalEnd = signalSet.end();
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = delayEmbed(*signalIter, k, t0, dt);
			++outputIter;
			++signalIter;
		}
	}

	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbedFuture(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer t0,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(t0, >=, 0);
		ENSURE_OP(dt, >=, 1);

		SignalPtr_Iterator signalIter = signalSet.begin();
		const SignalPtr_Iterator signalEnd = signalSet.end();
		OutputIterator outputIter = outputBegin;
		while(signalIter != signalEnd)
		{
			*outputIter = Tim::delayEmbedFuture(*signalIter, k, t0, dt);
			++outputIter;
			++signalIter;
		}
	}

}

#endif
