// Description: Mutual information between marginals of a normal distribution
// Documentation: mutual_information_analytic.txt

#ifndef TIM_MUTUAL_INFORMATION_NORMAL_H
#define TIM_MUTUAL_INFORMATION_NORMAL_H

#include "tim/core/mytypes.h"

#include <pastel/sys/real/real_concept.h>
#include <cmath>

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
		const NoDeduction<Real>& marginalCovarianceDeterminantProduct,
		const NoDeduction<Real>& jointCovarianceDeterminant)
	{
		/*
		The differential entropy of a multivariate gaussian
		with covariance C is given by:

		H(x) = 0.5 log((2 pi e)^d |C|)

		This is used to derive mutual information (total correlation)
		between x = (x_1, ..., x_m) and x_i's. Here x_i in R^(d_i),
		x in R^d, and sum_i d_i = d.
		
		MI = sum_i[H(x_i)] - H(x)
		= 0.5 sum_i[log((2 pi e)^d_i |C_i|)] - 0.5 log((2 pi e)^d |C|)
		= 0.5 sum_i log(|C_i|) + 0.5 sum_i[log((2 pi e)^d_i)] - 
		0.5 log((2 pi e)^d) - 0.5 log(|C|)
		= 0.5 sum_i log(|C_i|) - 0.5 log(|C|)
		= 0.5 log((|C_1| * ... * |C_m|) / |C|)
		*/
		
		return (Real)0.5 * std::log(
			marginalCovarianceDeterminantProduct /
			jointCovarianceDeterminant);
	}

}

#endif
