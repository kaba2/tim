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

	TIM Signal delayEmbed(
		const Signal& signal,
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
	*/
	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbed(
		Signal_Input inputSet,
		Signal_Output output,
		integer k,
		integer dt = 1);

	//! Returns the future of a signal under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	TIM Signal delayEmbedFuture(
		const Signal& signal,
		integer dt = 1);

	//! Computes the futures of signals under a given delay-embedding.
	/*!
	Preconditions:
	dt >= 1

	dt:
	Embedding delay.
	*/
	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbedFuture(
		Signal_Input inputSet,
		Signal_Output output,
		integer dt = 1);

}

#include "tim/core/delay_embed.hpp"

#endif
