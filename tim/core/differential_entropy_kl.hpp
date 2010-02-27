#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_KL_HPP

#include "tim/core/differential_entropy_kl.h"
#include "tim/core/generic_entropy.h"

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
			integer acceptedSamples, 
			integer estimateSamples) const
		{
			estimate *= (real)dimension / acceptedSamples;
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
		typename Real_OutputIterator,
		typename NormBijection>
	integer temporalDifferentialEntropyKl(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		const KlDifferential_EntropyAlgorithm<NormBijection>
			entropyAlgorithm(normBijection);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			result,
			kNearest);
	}

	template <
		typename SignalPtr_Iterator, 
		typename Real_OutputIterator>
	integer temporalDifferentialEntropyKl(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		return Tim::temporalDifferentialEntropyKl(
			signalSet,
			timeWindowRadius,
			result,
			kNearest,
			Default_NormBijection());
	}

	template <
		typename Real_OutputIterator,
		typename NormBijection>
	integer temporalDifferentialEntropyKl(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const NormBijection& normBijection)
	{
		return Tim::temporalDifferentialEntropyKl(
			forwardRange(constantIterator(signal)),
			timeWindowRadius,
			result,
			kNearest,
			normBijection);
	}

	template <typename Real_OutputIterator>
	integer temporalDifferentialEntropyKl(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		return Tim::temporalDifferentialEntropyKl(
			signal,
			timeWindowRadius,
			result,
			kNearest, 
			Default_NormBijection());
	}

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropyKl(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
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
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		integer kNearest)
	{
		return Tim::differentialEntropyKl(
			signalSet,
			kNearest,
			Default_NormBijection());
	}

	template <typename NormBijection>
	real differentialEntropyKl(
		const SignalPtr& signal,
		integer kNearest,
		const NormBijection& normBijection)
	{
		return Tim::differentialEntropyKl(
			forwardRange(constantIterator(signal)),
			kNearest,
			normBijection);
	}

}

#endif
