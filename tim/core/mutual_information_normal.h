// Description: Mutual information between marginals of a normal distribution
// Documentation: mutual_information_analytic.txt

#ifndef TIM_MUTUAL_INFORMATION_NORMAL_H
#define TIM_MUTUAL_INFORMATION_NORMAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real_concept.h>

namespace Tim
{

	//! Mutual information between marginals of a normal distribution.
	/*!
	Let X be a normally distributed random variable consisting of marginal
	random variables (X_1, ..., X_m). Then the mutual information
	(total correlation) between X_i is given by:
	
	I(X_1, ..., X_m) = 0.5 log((|cov(X_1)| ... |cov(X_m)| / |cov(X)))

	marginalCovarianceDeterminantProduct:
	|cov(X_1)| ... |cov(X_m)|

	jointCovarianceDeterminant:
	|cov(X)|
	*/
	template <typename Real>
	Real mutualInformationNormal(
		const PASTEL_NO_DEDUCTION(Real)& marginalCovarianceDeterminantProduct,
		const PASTEL_NO_DEDUCTION(Real)& jointCovarianceDeterminant);

}

#include "tim/core/mutual_information_normal.hpp"

#endif
