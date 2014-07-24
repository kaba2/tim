#ifndef TIM_TSALLIS_ENTROPY_UNIFORM_HPP
#define TIM_TSALLIS_ENTROPY_UNIFORM_HPP

#include "tim/core/tsallis_entropy_uniform.h"
#include "tim/core/differential_entropy_uniform.h"

namespace Tim
{

	template <typename Real>
	Real uniformTsallisEntropy(
		const PASTEL_NO_DEDUCTION(Real)& q,
		const PASTEL_NO_DEDUCTION(Real)& supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		if (q == 1)
		{
			return uniformDifferentialEntropy<Real>(supportVolume);
		}

		// I = 1 / m(S)^(q - 1)

		Real I = inverse(std::pow(supportVolume, q - 1));

		// H_q(X) = (1 - I) / (q - 1)

		return (1 - I) / (q - 1);
	}

}

#endif
