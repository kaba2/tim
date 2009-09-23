// Description: Differential entropy estimation
// Detail: Stowell-Plumbley recursive partition estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_SP_H
#define TIM_DIFFERENTIAL_ENTROPY_SP_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/iteratorrange.h>

namespace Tim
{

	template <typename SignalPtr_Iterator>
	real differentialEntropySp(
		const ForwardRange<SignalPtr_Iterator>& signalSet);	

}

#include "tim/core/differential_entropy_sp.hpp"

#endif
