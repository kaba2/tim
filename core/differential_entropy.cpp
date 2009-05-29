#include "tim/core/differential_entropy.h"

namespace Tim
{

	TIMCORE real gaussianDifferentialEntropy(
		integer dimension, real variance)
	{
		PENSURE1(dimension > 0, dimension);
		PENSURE1(variance >= 0, variance);

		return ((real)dimension / 2) * std::log(
			variance * 2 * constantPi<real>() * constantNeper<real>());
	}

	TIMCORE real uniformDifferentialEntropy(real supportVolume)
	{
		PENSURE1(supportVolume >= 0, supportVolume);

		return std::log(supportVolume);
	}

}
