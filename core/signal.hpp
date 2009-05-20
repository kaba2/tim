#ifndef TIM_SIGNAL_HPP
#define TIM_SIGNAL_HPP

#include "tim/core/signal.h"

namespace Tim
{

	inline SignalPtr newSignal(integer dimension, integer size)
	{
		return SignalPtr(new Signal(dimension, size));
	}

}

#endif
