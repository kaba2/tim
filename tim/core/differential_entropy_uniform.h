// Description: Differential entropy of a uniform distribution

#ifndef TIM_DIFFERENTIAL_ENTROPY_UNIFORM_H
#define TIM_DIFFERENTIAL_ENTROPY_UNIFORM_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real_concept.h>

namespace Tim
{

	//! Differential entropy of a uniform distribution.
	/*!
	Preconditions:
	supportVolume > 0

	supportVolume:
	The measure m of the support of the probability density
	function p: m({x in R^n: p(x) != 0})
	*/
	template <typename Real>
	Real differentialEntropyUniform(
		const PASTEL_NO_DEDUCTION(Real)& supportVolume);

}

#include "tim/core/differential_entropy_uniform.hpp"

#endif
