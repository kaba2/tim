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
				, jointEstimate(digamma<dreal>(kNearest_) - digamma<dreal>(n_))
				, baseEstimate(jointEstimate - std::log((dreal)kNearest_))
			{
			}

			dreal localJointEstimate() const
			{
				return jointEstimate;
			}

			dreal localMarginalEstimate(integer k) const
			{
				//return baseEstimate + std::log((dreal)k);
				return jointEstimate + std::log(
					1 + ((dreal)(k - kNearest) / n) * std::exp(-jointEstimate));
			}
	
		private:
			integer kNearest;
			integer n;
			dreal jointEstimate;
			dreal baseEstimate;
		};
	};

}

#endif
