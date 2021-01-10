// Description: Tsallis entropy estimation
// Detail: Leonenko-Pronzato-Savani nearest-neighbor estimator

#ifndef TIM_TSALLIS_ENTROPY_LPS_H
#define TIM_TSALLIS_ENTROPY_LPS_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include "pastel/math/normbijection/euclidean_normbijection.h"

#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"
#include "tim/core/differential_entropy_kl.h"

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
	inline TIM dreal tsallisDecideK(dreal q, integer kNearestSuggestion)
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

	class LpsTsallis_EntropyAlgorithm
	{
	public:
		// This algorithm computes the 
		// Leonenko-Pronzato-Savani estimator 
		// for Tsallis entropy.

		typedef Euclidean_Norm<dreal> Norm;

		LpsTsallis_EntropyAlgorithm(
			integer dimension,
			integer kNearest,
			dreal q)
			: distancePower_(dimension * (1 - q))
			, gammaRatio_(
				gamma<dreal>(kNearest) / 
				gamma<dreal>(kNearest + 1 - q))
			, finalFactor_(inverse(q - 1))
			, q_(q)
		{
		}

		const Norm& norm() const
		{
			return norm_;
		}
		
		template <Distance_Concept Distance>
		dreal sumTerm(Distance distance) const
		{
			// Let
			// C_k = (Gamma(k) / Gamma(k + 1 - q))^(1 / (1 - q))
			// V_m = Volume of unit ball in R^m
			// M = The total number of points participating
			//     in the estimate.
			// N = The number of non-zero d_i.

			// I
			// = (1 / N) sum_{i = 1}^M ((M - 1) C_k V_m (d_i)^m)^(1 - q)
			// = ((M - 1) C_k V_m)^(1 - q) (1 / N) sum_{i = 1}^M d_i^(m * (1 - q))
			// = F (1 / N) sum_{i = 1}^M d_i^(m * (1 - q))
			//
			// where
			//
			// F = ((M - 1) C_k V_m)^(1 - q)

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
			// F = ((M - 1) C_k V_m)^(1 - q)
			//   = ((M - 1) V_m)^(1 - q) (gamma(k) / gamma(k + 1 - q))
			//
			// Finally, Tsallis entropy estimator (for q != 1) is given by:
			// H_q(X) = (1 - I) / (q - 1)
			
			dreal F = 
				std::pow(
				(estimateSamples - 1) * 
				std::exp(lnVolumeUnitSphere(norm_, dimension)), 1 - q_) *
				gammaRatio_;

			const dreal I = F * estimate;

			dreal tsallisEntropy = 
				(1 - I) * finalFactor_;

			return tsallisEntropy;
		}

	private:
		Norm norm_;
		dreal distancePower_;
		dreal gammaRatio_;
		dreal finalFactor_;
		dreal q_;
	};

	// Temporal Tsallis entropy
	// ------------------------

	//! Computes temporal Tsallis entropy of a signal.
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

	q:
	The exponent in the definition of Tsallis entropy.
	If q == 1, the result of temporalDifferentialEntropyKl() 
	is returned instead.
	If q < 1, the results have huge errors: you
	should not use this estimator for those values.

	kNearestSuggestion:
	A suggestion for the k:th nearest neighbor that should be
	used for estimation. The k can't be set	freely because the 
	estimation algorithm is only defined for k > q - 1. 
	Value zero means an accurate (q-dependent) default is used.
	The actual k that is used is given by tsallisDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.
	*/
	template <
		ranges::forward_range Signal_Range, 
		ranges::forward_range Filter_Range>
	SignalData temporalTsallisEntropyLps(
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

		if (ranges::empty(signalSet))
		{
			return SignalData();
		}

		if (q == 1)
		{
			// Tsallis entropy is not defined for q = 1.
			// However, the limit of Tsallis entropy as q
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

		integer kNearest = tsallisDecideK(q, kNearestSuggestion);
		integer dimension = std::begin(signalSet)->dimension();
		
		LpsTsallis_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	//! Computes temporal Tsallis entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalTsallisEntropyLps(
		signalSet, timeWindowRadius,
		q, kNearestSuggestion,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
	template <ranges::forward_range Signal_Range>
	SignalData temporalTsallisEntropyLps(
		const Signal_Range& signalSet,
		integer timeWindowRadius,
		dreal q = 2,
		integer kNearestSuggestion = 0)
	{
		return Tim::temporalTsallisEntropyLps(
			signalSet, timeWindowRadius, 
			q, kNearestSuggestion,
			constantRange((dreal)1, 1));
	}

	// Tsallis entropy
	// ---------------

	//! Computes Tsallis entropy of a signal.
	/*!
	Preconditions:
	signalSet contains Signal's.
	q > 0
	kNearestSuggestion >= 0

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
	The actual k that is used is given by tsallisDecideK().
	For accurate results one should choose 
	kNearestSuggestion >= 2 * ceil(q) - 1.

	Returns:
	A Tsallis entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/
	template <ranges::forward_range Signal_Range>
	dreal tsallisEntropyLps(
		const Signal_Range& signalSet,
		dreal q = 2,
		integer kNearestSuggestion = 0)
	{
		ENSURE_OP(q, >, 0);
		ENSURE_OP(kNearestSuggestion, >=, 0);

		if (ranges::empty(signalSet))
		{
			return 0;
		}

		if (q == 1)
		{
			// Tsallis entropy is not defined for q = 1.
			// However, the limit of Tsallis entropy as q
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

		integer kNearest = tsallisDecideK(q, kNearestSuggestion);
		integer dimension = std::begin(signalSet)->dimension();
		
		LpsTsallis_EntropyAlgorithm entropyAlgorithm(
			dimension, kNearest, q);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			kNearest);
	}

}

#endif
