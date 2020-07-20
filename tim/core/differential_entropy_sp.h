// Description: Differential entropy estimation
// Detail: Stowell-Plumbley recursive partition estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_SP_H
#define TIM_DIFFERENTIAL_ENTROPY_SP_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/range.h>

namespace Tim
{

	template <ranges::forward_range Signal_Range>
	dreal differentialEntropySp(
		const Signal_Range& signalSet);	

}

#include "tim/core/differential_entropy_sp.hpp"

#endif
