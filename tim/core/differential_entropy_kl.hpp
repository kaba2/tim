#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_KL_HPP

#include "tim/core/differential_entropy_kl.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"

namespace Tim
{

	template <typename Used_NormBijection>
	class KlDifferential_EntropyAlgorithm
	{
	public:
		// This algorithm computes the Kozachenko-Leonenko
		// estimator for differential entropy.

		typedef Used_NormBijection NormBijection;

		KlDifferential_EntropyAlgorithm()
			: normBijection_()
		{
		}
		
		explicit KlDifferential_EntropyAlgorithm(
			const NormBijection& normBijection)
			: normBijection_(normBijection)
		{
		}

		const NormBijection& normBijection() const
		{
			return normBijection_;
		}
		
		real sumTerm(real distance) const
		{
			return normBijection_.toLnNorm(distance);
		}

		real finishEstimate(
			real estimate, 
			integer dimension, 
			integer kNearest, 
			integer estimateSamples) const
		{
			estimate *= (real)dimension;
			estimate -= digamma<real>(kNearest);
			estimate += digamma<real>(estimateSamples);
			estimate += normBijection_.lnVolumeUnitSphere(dimension);
			
			return estimate;
		}

	private:
		NormBijection normBijection_;
	};

	// Temporal differential entropy
	// -----------------------------

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection,
		typename Real_Filter_Iterator>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection,
		const ForwardIterator_Range<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		const KlDifferential_EntropyAlgorithm<NormBijection>
			entropyAlgorithm(normBijection);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection)
	{
		return Tim::temporalDifferentialEntropyKl(
			signalSet,
			timeWindowRadius,
			kNearest,
			normBijection,
			constantRange((real)1, 1));
	}

	template <typename SignalPtr_Iterator>
	SignalPtr temporalDifferentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		integer kNearest)
	{
		return Tim::temporalDifferentialEntropyKl(
			signalSet,
			timeWindowRadius,
			kNearest,
			Default_NormBijection());
	}

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(kNearest, >, 0);

		const KlDifferential_EntropyAlgorithm<NormBijection>
			entropyAlgorithm(normBijection);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			kNearest);
	}

	template <typename SignalPtr_Iterator>
	real differentialEntropyKl(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
		integer kNearest)
	{
		return Tim::differentialEntropyKl(
			signalSet,
			kNearest,
			Default_NormBijection());
	}

}

#endif
