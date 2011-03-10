// Description: Delay embedding

#ifndef TIM_DELAY_EMBED_H
#define TIM_DELAY_EMBED_H

#include "tim/core/signal.h"

#include <pastel/sys/iterator_range.h>

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
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer dt = 1);

	//! Time shift to produce the future of a delay-embedded signal.
	/*!
	Preconditions:
	k > 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.

	Assume a signal x : R -> R^n was embedded to a signal X : R -> R^(kn)
	using the function 'delayEmbed' with an embedding factor k. 
	Then the future of x is defined as the last R^n component of the 
	otherwise same embedding but with an embedding factor (k + 1).

	The future of x can be obtained simply by shifting x by
	a proper 'tDelta' amount of samples: 

	xFuture[t] = x[t + tDelta]
	
	This function computes 'tDelta' from the parameters that were used 
	for the delay embedding.
	*/

	TIM integer delayEmbedFutureShift(
		integer k, 
		integer dt = 1);

	//! Returns the future of a signal under a given delay-embedding.
	/*!
	Preconditions:
	k > 0
	t0 >= 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.
	*/
	TIM SignalPtr delayEmbedFuture(
		const SignalPtr& signal,
		integer k,
		integer dt = 1);

	//! Computes the futures of signals under a given delay-embedding.
	/*!
	Preconditions:
	k > 0
	dt >= 1

	k:
	Embedding factor.

	dt:
	Embedding delay.
	*/
	template <typename SignalPtr_Iterator, typename OutputIterator>
	void delayEmbedFuture(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		const OutputIterator& outputBegin,
		integer k,
		integer dt = 1);

}

#include "tim/core/delay_embed.hpp"

#endif
