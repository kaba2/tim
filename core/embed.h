// Description: Tools for delay embedding

#ifndef TIM_EMBED_H
#define TIM_EMBED_H

#include "tim/core/signal.h"

#include <vector>

namespace Tim
{

	//! Delay-embeds a signal into a higher dimensional space.
	/*!
	Preconditions:
	k > 0
	shift >= 0
	step >= 1

	Given is a signal S : Z -> R^n.
	Form a signal R : Z -> R^n : 
	R(t) = (S(t0 + t), S(t0 + t + dt), ..., S(t0 + t + dt * (k - 1)).
	
	Then R is the delay-embedding and
	t0 is the embedding shift
	dt is the embedding delay
	k is the 'embedding factor'
	d = k n is the embedding dimension
	*/

	TIMCORE SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer shift = 0,
		integer step = 1);

	//! Performs delay embedding for a set of signals.
	/*!
	Preconditions:
	k > 0
	shift >= 0
	step >= 1
	Signal_ForwardIterator dereferences to a SignalPtr.
	OutputIterator dereferences to a SignalPtr&.
	*/
	template <typename Signal_ForwardIterator, typename OutputIterator>
	void delayEmbed(
		const Signal_ForwardIterator& signalBegin,
		const Signal_ForwardIterator& signalEnd,
		const OutputIterator& outputBegin,
		integer k,
		integer shift = 0,
		integer step = 1);

	//! Time shift to produce the future of a delay-embedded signal.
	/*!
	Preconditions:
	k > 0
	shift >= 0
	step >= 1

	Assume that a signal x : R -> R^n was embedded to a signal X : R -> R^(kn)
	using the function 'delayEmbed' with an embedding factor k. 
	Then the future of x is defined as the last R^n component of the 
	otherwise same embedding but with an embedding factor (k + 1).

	The future of x can be obtained simply by shifting x by
	a proper 'tDelta' amount of samples: 

	xFuture[t] = x[t + tDelta]
	
	This function computes 'tDelta' from the parameters that were used 
	for the delay embedding.
	*/

	TIMCORE integer delayEmbedFutureShift(
		integer k, 
		integer shift = 0, 
		integer step = 1);


	TIMCORE SignalPtr delayEmbedFuture(
		const SignalPtr& signal,
		integer k,
		integer shift = 0,
		integer step = 1);

	template <typename Signal_ForwardIterator, typename OutputIterator>
	void delayEmbedFuture(
		const Signal_ForwardIterator& signalBegin,
		const Signal_ForwardIterator& signalEnd,
		const OutputIterator& outputBegin,
		integer k,
		integer shift = 0,
		integer step = 1);

}

#include "tim/core/embed.hpp"

#endif
