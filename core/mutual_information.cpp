#include "tim/core/mutual_information.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/histogram.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/matrix_tools.h>
#include <pastel/math/cholesky_decomposition_tools.h>

#include <pastel/gfx/pcx.h>
#include <pastel/gfx/image_tools.h>

namespace Tim
{

	TIMCORE real correlatedGaussianMutualInformation(
		real marginalCovarianceDeterminantProduct,
		real jointCovarianceDeterminant)
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
		
		return 0.5 * std::log(
			marginalCovarianceDeterminantProduct /
			jointCovarianceDeterminant);
	}

	template <typename ConstIterator>
	real entropy(
		const ConstIterator& begin,
		const ConstIterator& end,
		real binSize)
	{
		real result = 0;
		ConstIterator iter = begin;
		while(iter != end)
		{
			const real value = *iter;
			if (value > 0)
			{
				result -= value * std::log(value) * binSize;
			}
			++iter;
		}
		
		return result;
	}

	TIMCORE void mutualInformationNaive(
		const SignalPtr& signal,
		integer bins,
		MatrixD& result)
	{
		/*
		We consider 'signal' as a set of 1d signals. Each such signal has
		a continuous pdf. We approximate the mutual information
		for two 1d signals (called pairwise mutual information) 
		from their samplings as follows:
		1) Compute the min-max range of a signal.
		2) Divide the min-max range uniformly into bins. Enumerate
		these bins as integers in [0, 'bins'[.
		3) Associate each real number with the bin it falls into, giving
		a piecewise-constant distribution.
		4) Repeat 1-3 for the other signal.
		5) Compute mutual information for the piecewise-constant distributions.
		*/

		const integer samples = signal->samples();
		const integer n = signal->dimension();

		result.setSize(n, n);

		VectorD minBound = min(transpose(signal->data()));
		VectorD maxBound = max(transpose(signal->data()));
		VectorD binExtent = (maxBound - minBound) / bins;
		
		// Extend the bin support by a half bin
		// to guarantee that all samples fall into
		// some bin. Strictly, this is not needed,
		// but we want to make this function
		// behave almost equivalent to the implementation
		// in EEGLAB for comparison purposes.

		minBound -= binExtent / 2;
		maxBound += binExtent / 2;
		binExtent = (maxBound - minBound) / bins;

		Array<2, real> marginalHistogram(bins, n);

		for (integer i = 0;i < n;++i)
		{
			computeHistogram(
				signal->data().rowBegin(i),
				signal->data().rowEnd(i),
				minBound[i],
				maxBound[i],
				bins,
				marginalHistogram.rowBegin(i));
		}

		Array<2, real> jointHistogram(bins, bins);

		for (integer i = 0;i < n;++i)
		{
			for (integer j = i + 1;j < n;++j)
			{
				computeJointHistogram(
					signal->data().rowBegin(i),
					signal->data().rowEnd(i),
					minBound[i],
					maxBound[i],
					signal->data().rowBegin(j),
					signal->data().rowEnd(j),
					minBound[j],
					maxBound[j],
					arrayView(jointHistogram));

				const real binVolume = binExtent[i] * binExtent[j];

				real mi = 0;
				for (integer y = 0;y < bins;++y)
				{
					const real yMass = marginalHistogram(y, j);

					if (yMass > 0)
					{
						for (integer x = 0;x < bins;++x)
						{
							// We choose to do multiplications and division
							// instead of subtracting logarithms.
							// This way we hope to avoid possible cancellation
							// problems.

							const real xMass = marginalHistogram(x, i);
							const real xyMass = jointHistogram(x, y);
							if (xMass > 0 && xyMass > 0)
							{
								mi += xyMass * std::log(xyMass / (xMass * yMass)) * binVolume;
							}
						}
					}
				}

				result(i, j) = mi;
				result(j, i) = mi;
			}

			result(i, i) = entropy(
				marginalHistogram.rowBegin(i),
				marginalHistogram.rowEnd(i),
				binExtent[i]);
		}
	}

}

