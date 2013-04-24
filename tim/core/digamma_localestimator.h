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
				, nDigamma(digamma<real>(n_))
			{
			}

			real localJointEstimate() const
			{
				return digamma<real>(kNearest) - nDigamma;
			}

			real localMarginalEstimate(integer k) const
			{
				return digamma<real>(k) - nDigamma;
			}

		private:
			integer kNearest;
			real nDigamma;
		};
	};

}

#endif
