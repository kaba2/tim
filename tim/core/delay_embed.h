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

	TIM SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer dt = 1);

	//! Performs delay embedding for a set of signals.
	/*!
	Preconditions:
	k > 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.

	SignalPtr_Iterator dereferences to a SignalPtr.
	OutputIterator dereferences to a SignalPtr&.
	*/
	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbed(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer dt = 1);

	//! Returns the future of a signal under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	TIM SignalPtr delayEmbedFuture(
		const SignalPtr& signal,
		integer dt = 1);

	//! Computes the futures of signals under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbedFuture(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer dt = 1);

}

#include "tim/core/delay_embed.hpp"

#endif
