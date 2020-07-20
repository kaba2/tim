// Description: Kraskov-style expectation-of-log estimator
// Documentation: localestimators.txt

#ifndef TIM_LOG_LOCALESTIMATOR_H
#define TIM_LOG_LOCALESTIMATOR_H

#include "tim/core/localestimator_concept.h"

#include <pastel/sys/math_functions.h>

namespace Tim
{

	class Log_LocalEstimator
	{
	public:
		class Instance
		{
		public:
			Instance(
				integer kNearest_,
				integer n_)
				: kNearest(kNearest_)
				, nLog(std::log((dreal)n_))
			{
			}

			dreal localJointEstimate() const
			{
				return std::log((dreal)kNearest) - nLog;
			}

			dreal localMarginalEstimate(integer k) const
			{
				return std::log((dreal)k) - nLog;
			}

		private:
			integer kNearest;
			dreal nLog;
		};
	};

}

#endif
