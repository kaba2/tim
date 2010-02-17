#ifndef TIM_TSALLIS_ENTROPY_LPS_HPP
#define TIM_TSALLIS_ENTROPY_LPS_HPP

#include "tim/core/tsallis_entropy_lps.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/differential_entropy_kl.h"

#include "pastel/math/euclidean_normbijection.h"

namespace Tim
{

	class LpsTsallis_EntropyAlgorithm
	{
	public:
		// This algorithm computes the 
		// Leonenko-Pronzato-Savani estimator 
		// for Tsallis entropy.

		typedef Euclidean_NormBijection<real> NormBijection;

		LpsTsallis_EntropyAlgorithm(
			integer dimension,
			integer kNearest,
			real q)
			: distancePower_(dimension / (1 - q))
			, gammaRatio_(
				gamma<real>(kNearest) / 
				gamma<real>(kNearest + 1 - q))
			, finalFactor_(inverse(1 - q))
			, q_(q)
		{
		}

		const NormBijection& normBijection() const
		{
			return normBijection_;
		}
		
		real sumTerm(real distance) const
		{
			// Let
			// C_k = (Gamma(k) / Gamma(k + 1 - q))^(1 / (1 - q))
			// V_m = Volume of unit ball in R^m
			// M = The total number of points participating
			//     in the estimate.
			// N = The number of non-zero d_i.

			// I
			// = (1 / N) sum_{i = 1}^M ((M - 1) C_k V_m (d_i)^m)^(1 - q)
			// = ((M - 1) C_k V_m)^(1 - q) (1 / N) sum_{i = 1}^M d_i^(m / (1 - q))
			// = F (1 / N) sum_{i = 1}^M d_i^(m / (1 - q))
			//
			// where
			//
			// F = ((M - 1) C_k V_m)^(1 - q)

			return std::pow(
				normBijection_.toNorm(distance),
				distancePower_);
		}

		real finishEstimate(
			real estimate, 
			integer dimension, 
			integer kNearest, 
			integer acceptedSamples, 
			integer estimateSamples) const
		{
			// F = ((M - 1) C_k V_m)^(1 - q)
			//   = ((M - 1) V_m)^(1 - q) (gamma(k) / gamma(k + 1 - q))
			//
			// Finally, Tsallis entropy estimator (for q != 1) is given by:
			// H_q(X) = (1 - I) / (q - 1)
			
			const real F = 
				std::pow(
				(estimateSamples - 1) * 
				std::exp(normBijection_.lnVolumeUnitSphere(dimension)), 1 - q_) *
				gammaRatio_;

			const real tsallisEntropy = 
				finalFactor_ * (1 - F * (estimate / acceptedSamples));

			return tsallisEntropy;
		}

	private:
		NormBijection normBijection_;
		real distancePower_;
		real gammaRatio_;
		real finalFactor_;
		real q_;
	};

	// Temporal Tsallis entropy
	// ----------------------

	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator>
	bool temporalTsallisEntropyLps(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real q,
		real maxRelativeError,
		integer kNearest)
	{
		ENSURE_OP(timeWindowRadius, >=, 1);
		ENSURE_OP(maxRelativeError, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		if (q == 1)
		{
			// Tsallis entropy is not defined for q = 1.
			// However, the limit of Tsallis entropy as q
			// approaches 1 is given by the Shannon
			// differential entropy. We will return
			// this value instead for convenience.

			return temporalDifferentialEntropyKl(
				signalSet,
				timeWindowRadius,
				result,
				maxRelativeError,
				kNearest);
		}

		if (signalSet.empty())
		{
			return true;
		}

		const integer dimension = signalSet.front()->dimension();
		
		const LpsTsallis_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			result,
			maxRelativeError,
			kNearest);
	}

	template <typename Real_OutputIterator>
	bool temporalTsallisEntropyLps(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real q,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::temporalTsallisEntropyLps(
			signal,
			timeWindowRadius,
			result,
			q,
			maxRelativeError,
			kNearest);
	}

	// Tsallis entropy
	// -------------

	template <typename SignalPtr_Iterator>
	real tsallisEntropyLps(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		real q,
		real maxRelativeError,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		if (q == 1)
		{
			// Tsallis entropy is not defined for q = 1.
			// However, the limit of Tsallis entropy as q
			// approaches 1 is given by the Shannon
			// differential entropy. We will return
			// this value instead for convenience.

			return differentialEntropyKl(
				signalSet,
				maxRelativeError,
				kNearest);
		}

		const integer dimension = signalSet.front()->dimension();
		
		const LpsTsallis_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			maxRelativeError,
			kNearest);
	}

}

#endif
