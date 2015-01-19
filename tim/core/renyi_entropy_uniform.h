// Description: Renyi entropy of a uniform distribution
// Documentation: renyi_entropy_analytic.txt

#ifndef TIM_RENYI_ENTROPY_UNIFORM_H
#define TIM_RENYI_ENTROPY_UNIFORM_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real/real_concept.h>

namespace Tim
{

	//! Renyi entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume > 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})

	Note:
	The Renyi entropy of a uniform distribution
	is independent of q.
	*/
	template <typename Real>
	Real uniformRenyiEntropy(
		const NoDeduction<Real>& supportVolume);

}

#include "tim/core/renyi_entropy_uniform.hpp"

#endif
