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

}

#endif
