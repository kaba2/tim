#include "tim/core/renyi_entropy_analytic.h"
#include "tim/core/differential_entropy_analytic.h"

namespace Tim
{

	TIM real gaussianRenyiEntropy(
		real q,
		integer dimension,
		real covarianceDeterminant)
	{
		if (q == 1)
		{
			return gaussianDifferentialEntropy(
				dimension,
				covarianceDeterminant);
		}

		// H_q^*(X) = (n / 2)[log(2 pi) - (log(q) / (1 - q))] + 
		//            (1 / 2)log(det(C))

		return ((real)dimension / 2) * 
			(std::log(2 * constantPi<real>()) - std::log(q) / (1 - q)) +
			std::log(covarianceDeterminant) / 2;
	}

	TIM real uniformRenyiEntropy(
		real supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		return std::log(supportVolume);
	}

}
