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

		/*
		Let X be a random variable in R^n 
		with a probability density function p.
		Let S subset R^n be the support set of p.
		The measure of S is m(S). Then because
		of uniformity p(x) = 1 / m(S).
		Now
		H(X) = -int_S p(x) ln(p(x)) dx
		= -int_S (1 / m(s)) ln(1 / m(S)) dx
		= -int_S (1 / m(s)) (ln(1) - ln(m(S))) dx
		= -int_S (1 / m(s)) (-ln(m(S))) dx
		= (ln(m(S)) / m(S)) int_S 1 dx
		= (ln(m(S)) / m(S)) m(S)
		= ln(m(S))
		*/

		return std::log(supportVolume);
	}

}
