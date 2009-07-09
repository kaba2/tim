#ifndef TIM_EMBED_H
#define TIM_EMBED_H

#include "tim/core/signal.h"

#include <vector>

namespace Tim
{

	//! Delay-embeds a signal to a higher dimensional space.
	/*!
	Preconditions:
	k > 0
	shift >= 0
	step >= 1

	Let 
	t0 in Z be the embedding shift ('shift')
	dt in Z be the embedding step ('step')
	k in Z be the embedding factor ('k')
	x be a signal N -> R^n ('signal')

	Then the delay embedding of x is given by:
	y(t) = (x(t0 + dt t), x(t0 + dt (t + 1)), ..., x(t0 + dt (t + k - 1)))
	where
	y is a signal N -> R^(k n) (not (R^k)^n: the samples from x are concatenated).
	The embedding dimension is thus given by kn.
	*/

	TIMCORE SignalPtr delayEmbed(
		const SignalPtr& signal,
		integer k,
		integer shift = 0,
		integer step = 1);

}

#endif
