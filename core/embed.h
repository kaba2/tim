#ifndef TIM_EMBED_H
#define TIM_EMBED_H

#include "tim/core/signal.h"

#include <vector>

namespace Tim
{

	//! Delay-embeds a signal to a high-dimensional space.
	/*!
	Preconditions:
	dimension > 0
	step >= 1

	Delay embedding constructs a vector set 
	{v_1, ..., v_m} in R^n by
	v_i[j] = data[(i * n + j) * step]
	*/

	TIMCORE SignalPtr delayEmbed(
		const std::vector<real>& data,
		integer dimension,
		integer step = 1);

}

#endif
