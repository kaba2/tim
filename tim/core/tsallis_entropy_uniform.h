// Description: Tsallis entropy of a uniform distribution
// Documentation: tsallis_entropy_analytic.txt

#ifndef TIM_TSALLIS_ENTROPY_UNIFORM_H
#define TIM_TSALLIS_ENTROPY_UNIFORM_H

#include "tim/core/mytypes.h"
#include "tim/core/differential_entropy_uniform.h"

#include <pastel/sys/real/real_concept.h>

namespace Tim
{

	//! Tsallis entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume >= 0

	q:
	The power in the definition of Tsallis entropy.
	If q == 1, differentialEntropyUniform() is
	returned instead.

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})
	*/
	template <typename Real>
	Real uniformTsallisEntropy(
		const NoDeduction<Real>& q,
		const NoDeduction<Real>& supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		if (q == 1)
		{
			return differentialEntropyUniform<Real>(supportVolume);
		}

		// I = 1 / m(S)^(q - 1)

		Real I = inverse(std::pow(supportVolume, q - 1));

		// H_q(X) = (1 - I) / (q - 1)

		return (1 - I) / (q - 1);
	}

}

#include "tim/core/tsallis_entropy_uniform.hpp"

#endif
