// Description: My expectation-of-log estimator
// Documentation: localestimators.txt

#ifndef TIM_DIGAMMADENSITY_LOCALESTIMATOR_H
#define TIM_DIGAMMADENSITY_LOCALESTIMATOR_H

#include "tim/core/localestimator_concept.h"

#include <pastel/sys/math_functions.h>

namespace Tim
{

	class DigammaDensity_LocalEstimator
	{
	public:
		class Instance
		{
		public:
			Instance(
				integer kNearest_,
				integer n_)
				: kNearest(kNearest_)
				, n(n_)
				, jointEstimate(digamma<real>(kNearest_) - digamma<real>(n_))
				, baseEstimate(jointEstimate - std::log((real)kNearest_))
			{
			}

			real localJointEstimate() const
			{
				return jointEstimate;
			}

			real localMarginalEstimate(integer k) const
			{
				//return baseEstimate + std::log((real)k);
				return jointEstimate + std::log(
					1 + ((real)(k - kNearest) / n) * std::exp(-jointEstimate));
			}
	
		private:
			integer kNearest;
			integer n;
			real jointEstimate;
			real baseEstimate;
		};
	};

}

#endif
