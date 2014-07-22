#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_KL_HPP

#include "tim/core/differential_entropy_kl.h"
#include "tim/core/generic_entropy.h"
#include "tim/core/generic_entropy_t.h"

namespace Tim
{

	template <typename NormBijection_>
	class KlDifferential_EntropyAlgorithm
	{
	public:
		// This algorithm computes the Kozachenko-Leonenko
		// estimator for differential entropy.

		using NormBijection = NormBijection_;

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

	template <
		typename SignalPtr_Range, 
		typename NormBijection,
		typename Real_Range>
	Signal temporalDifferentialEntropyKl(
		const SignalPtr_Range& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const NormBijection& normBijection,
		const Real_Range& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		KlDifferential_EntropyAlgorithm<NormBijection>
			entropyAlgorithm(normBijection);

		return temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			filter);
	}

	template <
		typename SignalPtr_Range, 
		typename NormBijection>
	real differentialEntropyKl(
		const SignalPtr_Range& signalSet,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(kNearest, >, 0);

		KlDifferential_EntropyAlgorithm<NormBijection>
			entropyAlgorithm(normBijection);

		return genericEntropy(
			signalSet,
			entropyAlgorithm,
			kNearest);
	}

}

#endif
