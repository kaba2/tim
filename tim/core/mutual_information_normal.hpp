#ifndef TIM_MUTUAL_INFORMATION_NORMAL_HPP
#define TIM_MUTUAL_INFORMATION_NORMAL_HPP

#include "tim/core/mutual_information_normal.h"

#include <cmath>

namespace Tim
{

	template <typename Real>
	Real mutualInformationNormal(
		const PASTEL_NO_DEDUCTION(Real)& marginalCovarianceDeterminantProduct,
		const PASTEL_NO_DEDUCTION(Real)& jointCovarianceDeterminant)
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
