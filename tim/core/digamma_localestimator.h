// Description: Kraskov's log-of-expectation estimator
// Documentation: localestimators.txt

#ifndef TIM_DIGAMMA_LOCALESTIMATOR_H
#define TIM_DIGAMMA_LOCALESTIMATOR_H

#include "tim/core/localestimator_concept.h"

#include <pastel/sys/math_functions.h>

namespace Tim
{

	class Digamma_LocalEstimator
	{
	public:
		class Instance
		{
		public:
			Instance(
				integer kNearest_,
				integer n_)
				: kNearest(kNearest_)
				, nDigamma(digamma<dreal>(n_))
			{
			}

			dreal localJointEstimate() const
			{
				return digamma<dreal>(kNearest) - nDigamma;
			}

			dreal localMarginalEstimate(integer k) const
			{
				return digamma<dreal>(k) - nDigamma;
			}

		private:
			integer kNearest;
			dreal nDigamma;
		};
	};

}

#endif
