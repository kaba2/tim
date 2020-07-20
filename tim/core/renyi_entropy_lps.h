// Description: Renyi entropy estimation
// Detail: Leonenko-Pronzato-Savani nearest-neighbor estimator

#ifndef TIM_RENYI_ENTROPY_LPS_H
#define TIM_RENYI_ENTROPY_LPS_H

#include "tim/core/signal.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"
#include "tim/core/differential_entropy_kl.h"

#include <pastel/sys/range.h>
#include "pastel/math/normbijection/euclidean_normbijection.h"

namespace Tim
{

	//! Returns the actual k:th neighbor to use given a suggestion.
	/*!
	Preconditions:
	q > 0
	kNearestSuggestion >= 0

	Returns:
	If kNearestSuggestion == 0, then k = 2 * ceil(q).
	If 0 < kNearestSuggestion < q - 1, then k = ceil(q - 1).
	If kNearestSuggestion == q - 1, then k = kNearestSuggestion + 1.
	Otherwise, k = kNearestSuggestion.

	The k in the Leonenko-Pronzato-Savani estimator can't be set 
	freely because the algorithm is not defined when k <= q - 1.
	This functions helps to decide a proper k.
	*/
	inline TIM dreal renyiDecideK(dreal q, integer kNearestSuggestion)
	{
		PENSURE_OP(q, >, 0);
		PENSURE_OP(kNearestSuggestion, >=, 0);

		integer kNearest = kNearestSuggestion;

		if (kNearestSuggestion == 0)
		{
			// We get to decide the k.

			kNearest = 2 * std::ceil(q);
		}
		else if (kNearestSuggestion <= q - 1)
		{
			// The algorithm is not defined
			// for such k. Find the smallest
			// k for which the algorithm is
			// defined.

			if (kNearestSuggestion < q - 1)
			{
				// 0 < kNearestSuggestion
				// thus
				// 0 < q - 1

				kNearest = std::ceil(q - 1);
			}
			else
			{
				kNearest = kNearestSuggestion + 1;
			}
		}

		return kNearest;
	}


	class LpsRenyi_EntropyAlgorithm
	{
	public:
		// This algorithm computes the 
		// Leonenko-Pronzato-Savani estimator 
		// for Renyi entropy.

		typedef Euclidean_Norm<dreal> Norm;

		LpsRenyi_EntropyAlgorithm(
			integer dimension,
			integer kNearest,
			dreal q)
			: distancePower_(dimension * (1 - q))
			, lnGammaDifference_(
				lnGamma<dreal>(kNearest) -
				lnGamma<dreal>(kNearest + 1 - q))
			, finalFactor_(inverse(1 - q))
			, q_(q)
		{
		}

		const Norm& norm() const
		{
			return norm_;
		}
		
		dreal sumTerm(Distance_Concept auto distance) const
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
				(dreal)distance,
				distancePower_);
		}

		dreal finishEstimate(
			dreal estimate, 
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

			const dreal F = (1 - q_) * (
				std::log((dreal)(estimateSamples - 1)) + 
				lnVolumeUnitSphere(norm_, dimension)) +
				lnGammaDifference_;
			
			dreal renyiEntropy = 
				finalFactor_ * (F + std::log(estimate));

			return renyiEntropy;
		}

	private:
		Norm norm_;
		dreal distancePower_;
		dreal lnGammaDifference_;
		dreal finalFactor_;
		dreal q_;
	};

	// Temporal Renyi entropy
	// -----------------------------

	//! Computes temporal Renyi entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearestSuggestion >= 0
	q > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	result:
	A dreal output iterator, denoting the start
	of the region where the sequence of temporal
	Renyi entropies are to be stored.

	q:
	The exponent in the definition of Renyi entropy.
	If q == 1, the result of temporalDifferentialEntropyKl() 
	is returned instead.
	If q < 1, the results have huge errors: you
	should not use this estimator for those values.

	kNearestSuggestion:
	A suggestion for the k:th nearest neighbor that should be
	used for estimation. The k can't be set	freely because the 
	estimation algorithm is only defined for k > q - 1. 
	Value zero means an accurate (q-dependent) default is used.
	The actual k that is used is given by renyiDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.

	Returns:
	The number of time instants that had an
	undefined estimate. If not all estimates
	were undefined, they were reconstructed from 
	the defined estimates using interpolation.
	*/
	template <
		ranges::forward_range Signal_Range, 
		ranges::forward_range Filter_Range>
	SignalData temporalRenyiEntropyLps(
		const Signal_Range& signalSet,
		integer timeWindowRadius,
		dreal q,
		integer kNearestSuggestion,
		const Filter_Range& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(q, >, 0);
		ENSURE_OP(kNearestSuggestion, >=, 0);
		ENSURE(odd(ranges::size(filter)));

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

		if (ranges::empty(signalSet))
		{
			return SignalData(0, 1);
		}

		integer kNearest = renyiDecideK(q, kNearestSuggestion);
		integer dimension = std::begin(signalSet)->dimension();
		
		LpsRenyi_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	//! Computes temporal Renyi entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalRenyiEntropyLps(
		signalSet, timeWindowRadius, result, 
		q, kNearestSuggestion,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
	template <ranges::forward_range Signal_Range>
	Signal temporalRenyiEntropyLps(
		const Signal_Range& signalSet,
		integer timeWindowRadius,
		dreal q = 2,
		integer kNearestSuggestion = 0)
	{
		return Tim::temporalRenyiEntropyLps(
			signalSet, timeWindowRadius,
			q, kNearestSuggestion,
			constantRange((dreal)1, 1));
	}

	// Renyi entropy
	// --------------------

	//! Computes Renyi entropy of a signal.
	/*!
	Preconditions:
	kNearestSuggestion > 0
	signalSet contains Signal's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	q:
	The exponent in the definition of Tsallis entropy.
	If q == 1, the result of differentialEntropyKl() 
	is returned instead.
	If q < 1, the results have huge errors: you
	should not use this estimator for those values.

	kNearestSuggestion:
	A suggestion for the k:th nearest neighbor that should be
	used for estimation. The k can't be set	freely because the 
	estimation algorithm is only defined for k > q - 1. 
	Value zero means an accurate (q-dependent) default is used.
	The actual k that is used is given by renyiDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.

	Returns:
	A Renyi entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/
	template <ranges::forward_range Signal_Range>
	dreal renyiEntropyLps(
		const Signal_Range& signalSet,
		dreal q = 2,
		integer kNearestSuggestion = 0)
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

		integer kNearest = renyiDecideK(q, kNearestSuggestion);
		integer dimension = std::begin(signalSet)->dimension();
		
		LpsRenyi_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			kNearest);
	}

}

#endif
