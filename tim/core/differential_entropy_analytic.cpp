#include "tim/core/differential_entropy_analytic.h"

namespace Tim
{

	TIM real gaussianDifferentialEntropy(
		integer dimension, real covarianceDeterminant)
	{
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(covarianceDeterminant, >=, 0);

		static const real ConstantFactor = std::log(
			2 * constantPi<real>()) + 1;

		return 0.5 * (std::log(std::abs(covarianceDeterminant)) + 
			dimension * ConstantFactor);
	}

	TIM real uniformDifferentialEntropy(real supportVolume)
	{
		PENSURE_OP(supportVolume, >, 0);

		// Let X be a random variable in R^n 
		// with a probability density function p.
		// Let S subset R^n be the support set of p.
		// The measure of S is m(S). Then because
		// of uniformity p(x) = 1 / m(S).
		// Now
		// H(X) = -int_S p(x) ln(p(x)) dx
		// = -int_S (1 / m(s)) ln(1 / m(S)) dx
		// = -int_S (1 / m(s)) (ln(1) - ln(m(S))) dx
		// = -int_S (1 / m(s)) (-ln(m(S))) dx
		// = (ln(m(S)) / m(S)) int_S 1 dx
		// = (ln(m(S)) / m(S)) m(S)
		// = ln(m(S))

		return std::log(supportVolume);
	}

	TIM real generalizedGaussianDifferentialEntropy(
		integer dimension, real shape, real scale)
	{
		PENSURE_OP(dimension, >, 0);
		PENSURE_OP(shape, >, 0);
		PENSURE_OP(scale, >, 0);

		// Let
		// a = scale
		// b = shape
		//
		// Then
		//
		// differential entropy 
		// = (1 / b) - log(b / (2a gamma(1 / b)))
		// = (1 / b) - (log(b / (2a)) - log(gamma(1 / b)))
		// = (1 / b) - log(b / (2a)) + log(gamma(1 / b))
		// = (1 / b) + log((2a) / b) + log(gamma(1 / b))
		//
		// The point of this derivation is to evaluate lnGamma 
		// instead of gamma for better numerical behaviour.

		real invShape = 1 / shape;


		return dimension * (invShape + std::log(2 * scale * invShape) +
			lnGamma<real>(invShape));
	}

}
