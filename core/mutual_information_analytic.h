// Description: Analytic solutions for mutual information

#ifndef TIM_MUTUAL_INFORMATION_ANALYTIC_H
#define TIM_MUTUAL_INFORMATION_ANALYTIC_H

#include "tim/core/mytypes.h"

namespace Tim
{

	//! Mutual information between marginals of a correlated gaussian.
	/*!
	Let X be a correlated gaussian random variable consisting of marginal
	random variables (X_1, ..., X_m). Then the mutual information
	(total correlation) between X_i is given by:
	
	I(X_1, ..., X_m) = 0.5 (|cov(X_1)| ... |cov(X_m)| / |cov(X))

	marginalCovarianceDeterminantProduct:
	|cov(X_1)| ... |cov(X_m)|

	jointCovarianceDeterminant:
	|cov(X)|
	*/

	TIMCORE real correlatedGaussianMutualInformation(
		real marginalCovarianceDeterminantProduct,
		real jointCovarianceDeterminant);

}

#endif
