#ifndef TIM_DELAY_EMBED_HPP
#define TIM_DELAY_EMBED_HPP

#include "tim/core/delay_embed.h"

namespace Tim
{

	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbed(
		Signal_Input inputSet,
		Signal_Output output,
		integer k,
		integer dt)
	{
		ENSURE_OP(k, >, 0);
		ENSURE_OP(dt, >=, 1);

		while(!inputSet.empty())
		{
			output(delayEmbed(inputSet.get(), k, dt));
			inputSet.pop();
		}
	}

	template <
		typename Signal_Input, 
		typename Signal_Output>
	void delayEmbedFuture(
		Signal_Input inputSet,
		Signal_Output output,
		integer dt)
	{
		ENSURE_OP(dt, >=, 1);

		while(!inputSet.empty())
		{
			output(Tim::delayEmbedFuture(inputSet.get(), dt));
			inputSet.pop();
		}
	}

}

#endif
