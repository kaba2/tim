// Description: Some analytic solutions of mutual information

#ifndef TIM_MUTUAL_INFORMATION_ANALYTIC_H
#define TIM_MUTUAL_INFORMATION_ANALYTIC_H

#include "tim/core/mytypes.h"

namespace Tim
{

	TIMCORE real correlatedGaussianMutualInformation(
		real marginalCovarianceDeterminantProduct,
		real jointCovarianceDeterminant);

}

#endif
