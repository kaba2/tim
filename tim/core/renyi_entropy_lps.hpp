#ifndef TIM_RENYI_ENTROPY_LPS_HPP
#define TIM_RENYI_ENTROPY_LPS_HPP

#include "tim/core/renyi_entropy_lps.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"
#include "tim/core/differential_entropy_kl.h"

#include "pastel/math/euclidean_normbijection.h"

namespace Tim
{

	class LpsRenyi_EntropyAlgorithm
	{
	public:
		// This algorithm computes the 
		// Leonenko-Pronzato-Savani estimator 
		// for Renyi entropy.

		typedef Euclidean_NormBijection<real> NormBijection;

		LpsRenyi_EntropyAlgorithm(
			integer dimension,
			integer kNearest,
			real q)
			: distancePower_(dimension * (1 - q))
			, lnGammaDifference_(
				lnGamma<real>(kNearest) -
				lnGamma<real>(kNearest + 1 - q))
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

			// log(I) 
			// = log[(1 / N) sum_{i = 1}^M ((M - 1) C_k V_m (d_i)^m)^(1 - q)]
			// = log[((M - 1) C_k V_m)^(1 - q) (1 / N) sum_{i = 1}^M d_i^(m(1 - q))]
			// = F + log[(1 / N) sum_{i = 1}^M d_i^(m(1 - q))]
			//
			// where
			//
			// F = log[((M - 1) C_k V_m)^(1 - q)]

			return std::pow(
				normBijection_.toNorm(distance),
				distancePower_);
		}

		real finishEstimate(
			real estimate, 
			integer dimension, 
			integer kNearest, 
			integer estimateSamples) const
		{
			// F = log[((M - 1) C_k V_m)^(1 - q)]
			// = (1 - q) log[(M - 1) C_k V_m]
			// = (1 - q) (log[M - 1] + log[V_m]) +
			//   (1 - q) log[C_k]
			// = (1 - q) (log[M - 1] + log[V_m]) +
			//   (lnGamma(k) - lnGamma(k + 1 - q))
			//
			// Finally, Renyi entropy estimator (for q != 1) is given by:
			// H_q(X) = (1 / (1 - q)) log(I)

			const real F = (1 - q_) * (
				std::log((real)(estimateSamples - 1)) + 
				normBijection_.lnVolumeUnitSphere(dimension)) +
				lnGammaDifference_;
			
			const real renyiEntropy = 
				finalFactor_ * (F + std::log(estimate));

			return renyiEntropy;
		}

	private:
		NormBijection normBijection_;
		real distancePower_;
		real lnGammaDifference_;
		real finalFactor_;
		real q_;
	};

	// Temporal Renyi entropy
	// ----------------------

	template <
		typename SignalPtr_Iterator, 
		typename Real_Filter_Iterator>
	SignalPtr temporalRenyiEntropyLps(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		real q,
		integer kNearestSuggestion,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(q, >, 0);
		ENSURE_OP(kNearestSuggestion, >=, 0);
		ENSURE(odd(filter.size()));

		if (q == 1)
		{
			// Renyi entropy is not defined for q = 1.
			// However, the limit of Renyi entropy as q
			// approaches 1 is given by the Shannon
			// differential entropy. We will return
			// this value instead for convenience.

			integer kNearest = kNearestSuggestion;
			if (kNearestSuggestion == 0)
			{
				kNearest = 1;
			}

			return temporalDifferentialEntropyKl(
				signalSet,
				timeWindowRadius,
				kNearest);
		}

		if (signalSet.empty())
		{
			return SignalPtr(new Signal(0, 1));
		}

		const integer kNearest = renyiDecideK(q, kNearestSuggestion);
		const integer dimension = signalSet.front()->dimension();
		
		const LpsRenyi_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	template <typename SignalPtr_Iterator>
	SignalPtr temporalRenyiEntropyLps(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		real q,
		integer kNearestSuggestion)
	{
		return Tim::temporalRenyiEntropyLps(
			signalSet, timeWindowRadius,
			q, kNearestSuggestion,
			constantRange((real)1, 1));
	}

	// Renyi entropy
	// -------------

	template <typename SignalPtr_Iterator>
	real renyiEntropyLps(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		real q,
		integer kNearestSuggestion)
	{
		ENSURE_OP(q, >, 0);
		ENSURE_OP(kNearestSuggestion, >=, 0);

		if (q == 1)
		{
			// Renyi entropy is not defined for q = 1.
			// However, the limit of Renyi entropy as q
			// approaches 1 is given by the Shannon
			// differential entropy. We will return
			// this value instead for convenience.

			integer kNearest = kNearestSuggestion;
			if (kNearestSuggestion == 0)
			{
				kNearest = 1;
			}

			return differentialEntropyKl(
				signalSet,
				kNearest);
		}

		const integer kNearest = renyiDecideK(q, kNearestSuggestion);
		const integer dimension = signalSet.front()->dimension();
		
		const LpsRenyi_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			kNearest);
	}

}

#endif
