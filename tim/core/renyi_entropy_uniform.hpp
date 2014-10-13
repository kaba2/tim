#ifndef TIM_RENYI_ENTROPY_UNIFORM_HPP
#define TIM_RENYI_ENTROPY_UNIFORM_HPP

#include "tim/core/renyi_entropy_uniform.h"

#include <cmath>

namespace Tim
{

	template <typename Real>
	Real uniformRenyiEntropy(
		const NoDeduction<Real>& supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		return std::log(supportVolume);
	}

}

#endif
