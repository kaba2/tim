#include "tim/core/tsallis_entropy_analytic.h"
#include "tim/core/differential_entropy_analytic.h"

namespace Tim
{

	TIM real gaussianTsallisEntropy(
		real q,
		integer dimension,
		real covarianceDeterminant)
	{
		PENSURE_OP(q, >, 0);
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(covarianceDeterminant, >=, 0);

		if (q == 1)
		{
			return gaussianDifferentialEntropy(
				dimension,
				covarianceDeterminant);
		}

		// A = 1 / ((2 pi)^(n / 2) det(C)^(1 / 2))
		
		real A = inverse(

			std::pow(2 * constantPi<real>(), (real)dimension / 2) *
			std::sqrt(covarianceDeterminant));

		// I = A^(q - 1) / q^(n / 2)

		real I = std::pow(A, q - 1) / 
			std::pow(q, (real)dimension / 2);

		// H_q(X) = (1 - I) / (q - 1)

		return (1 - I) / (q - 1);
	}

	TIM real uniformTsallisEntropy(
		real q,
		real supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		if (q == 1)
		{
			return uniformDifferentialEntropy(supportVolume);
		}

		// I = 1 / m(S)^(q - 1)

		real I = inverse(std::pow(supportVolume, q - 1));

		// H_q(X) = (1 - I) / (q - 1)

		return (1 - I) / (q - 1);
	}

}
