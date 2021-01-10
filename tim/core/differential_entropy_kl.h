// Description: Differential entropy estimation
// Detail: Kozachenko-Leonenko nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_H
#define TIM_DIFFERENTIAL_ENTROPY_KL_H

#include "tim/core/signal.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"

#include <pastel/sys/range.h>

#include <pastel/math/normbijection/normbijection_concept.h>

namespace Tim
{

	template <typename Norm_>
	class KlDifferential_EntropyAlgorithm
	{
	public:
		// This algorithm computes the Kozachenko-Leonenko
		// estimator for differential entropy.

		using Norm = Norm_;

		KlDifferential_EntropyAlgorithm()
			: norm_()
		{
		}
		
		explicit KlDifferential_EntropyAlgorithm(
			const Norm& norm)
			: norm_(norm)
		{
		}

		const Norm& norm() const
		{
			return norm_;
		}
		
		template <Distance_Concept Distance>
		dreal sumTerm(Distance distance) const
		{
			return std::log((dreal)distance);
		}

		dreal finishEstimate(
			dreal estimate, 
			integer dimension, 
			integer kNearest, 
			integer estimateSamples) const
		{
			estimate *= (dreal)dimension;
			estimate -= digamma<dreal>(kNearest);
			estimate += digamma<dreal>(estimateSamples);
			estimate += lnVolumeUnitSphere(norm_, dimension);
			
			return estimate;
		}

	private:
		Norm norm_;
	};

}

namespace Tim
{

	//! Temporal differential entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	norm:
	The norm to use.
	*/
	template <
		ranges::forward_range Signal_Range, 
		typename Norm = Default_Norm,
		typename Real_Range = decltype(constantRange((dreal)1, 1))>
	SignalData temporalDifferentialEntropyKl(
		const Signal_Range& signalSet,
		integer timeWindowRadius,
		integer kNearest = 1,
		const Norm& norm = Norm(),
		const Real_Range& filter = constantRange((dreal)1, 1))
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		KlDifferential_EntropyAlgorithm<Norm>
			entropyAlgorithm(norm);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	//! Differential entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0
	signalSet contains Signal's.

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate differential entropy.

	norm:
	The norm to use.

	Returns:
	A differential entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/
	template <
		ranges::forward_range Signal_Range, 
		typename Norm = Default_Norm>
	dreal differentialEntropyKl(
		const Signal_Range& signalSet,
		integer kNearest = 1,
		const Norm& norm = Norm())
	{
		ENSURE_OP(kNearest, >, 0);

		KlDifferential_EntropyAlgorithm<Norm> entropyAlgorithm(norm);
		return genericEntropy(signalSet, entropyAlgorithm, kNearest);
	}

}

#endif
