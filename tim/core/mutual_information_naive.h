// Description: Mutual information estimation using naive algorithms
// Detail: Includes computation from differential entropies and from a binning.

#ifndef TIM_MUTUAL_INFORMATION_NAIVE_H
#define TIM_MUTUAL_INFORMATION_NAIVE_H

#include "tim/core/signal.h"

#include <pastel/sys/histogram.h>

#include <pastel/math/matrix/matrix.h>
#include <pastel/sys/view/arrayview.h>

namespace Tim
{

	//! Computes mutual information from entropies.
	/*!
	Preconditions:
	kNearest > 0

	signalSet:
	A set of signals between which the mutual information
	is computed.

	kNearest:
	The k:th neighbor to use in the estimation of the
	differential entropies.

	norm:
	The norm to use to do the estimations.
	See 'pastel/math/normbijection/normbijections.h'.

	This technique is not recommended because it tends
	to give large errors in estimation. The intent of
	this function is to demonstrate the non-applicability
	of the technique.
	*/
	template <
		ranges::forward_range Signal_Range,
		typename Norm,
		typename Real_OutputIterator>
	dreal mutualInformationFromEntropy(
		const Signal_Range& signalSet,
		integer timeWindowRadius,
		integer kNearest,
		const Norm& norm,
		Real_OutputIterator result)
	{
		ENSURE_OP(kNearest, >, 0);

		if (ranges::empty(signalSet))
		{
			return 0;
		}

		auto iter = signalSet.begin();
		auto iterEnd = signalSet.end();

		dreal estimate = 0;
		while(iter != iterEnd)
		{
			const Signal& signal = *iter;

			estimate += differentialEntropyKl(
				signal,
				kNearest,
				norm);

			++iter;
		}

		const Signal& jointSignal = merge(signalSet);
		estimate -= differentialEntropyKl(
			jointSignal,
			kNearest,
			norm);

		return estimate;
	}

}

namespace Tim
{

	namespace
	{

		template <typename ConstIterator>
		dreal entropy(
			const ConstIterator& begin,
			const ConstIterator& end,
			dreal binSize)
		{
			dreal result = 0;
			ConstIterator iter = begin;
			while(iter != end)
			{
				const dreal value = *iter;
				if (value > 0)
				{
					result -= value * std::log(value) * binSize;
				}
				++iter;
			}
			
			return result;
		}

	}

	//! Computes pairwise 1d mutual information by binning.
	/*!
	Preconditions:
	bins > 0

	signal:
	Signal that contains n

	bins:
	Number of bins to use to estimate a 1d probability
	distribution function (pdf). For the joint pdf, a
	2d array of extents bins x bins is used.

	result (output):
	The element (i, j) contains the mutual information
	between the i:th and j:th 1d marginal signals of
	the 'signal'.

	The approximation of probability distribution functions
	using binning does not generalize practically to higher
	dimensions than 2, because of the exponential explosion
	of the number of needed bins. Therefore, using this
	technique, the mutual information can only be computed
	between two 1d signals. Given a multi-dimensional 
	signal, this function computes pairwise mutual
	information between the 1d marginal signals of the
	given signal.

	This technique is not recommended because it tends
	to give large errors in estimation. The intent of
	this function is to demonstrate the non-applicability
	of the technique.
	*/
	inline TIM void mutualInformationFromBinning(
		const Signal& signal,
		integer bins,
		const MatrixView<dreal>& result)
	{
		ENSURE_OP(bins, >, 0);

		/*
		We consider 'signal' as a set of 1d signals. Each such signal has
		a continuous pdf. We approximate the mutual information
		for two 1d signals (called pairwise mutual information) 
		from their samplings as follows:
		1) Compute the min-max range of a signal.
		2) Divide the min-max range uniformly into bins. Enumerate
		these bins as integers in [0, 'bins'[.
		3) Associate each dreal number with the bin it falls into, giving
		a piecewise-constant distribution.
		4) Repeat 1-3 for the other signal.
		5) Compute mutual information for the piecewise-constant distributions.
		*/

		integer n = signal.dimension();

		ENSURE_OP(result.rows(), ==, n);
		ENSURE_OP(result.cols(), ==, n);

		VectorD minBound = asVector(asMatrix(signal.data()).colwise().minCoeff());
		VectorD maxBound = asVector(asMatrix(signal.data()).colwise().maxCoeff());
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

		Array<dreal> marginalHistogram(Vector2i(bins, n));

		for (integer i = 0;i < n;++i)
		{
			computeHistogram(
				signal.data().columnRange(i),
				minBound[i],
				maxBound[i],
				bins,
				marginalHistogram.rowBegin(i));
		}

		Array<dreal> jointHistogram(Vector2i(bins, bins));

		for (integer i = 0;i < n;++i)
		{
			for (integer j = i + 1;j < n;++j)
			{
				computeJointHistogram(
					std::begin(signal.data().columnRange(i)),
					std::end(signal.data().columnRange(i)),
					minBound[i],
					maxBound[i],
					std::begin(signal.data().columnRange(j)),
					std::end(signal.data().columnRange(j)),
					minBound[j],
					maxBound[j],
					arrayView(jointHistogram));


				const dreal binVolume = binExtent[i] * binExtent[j];

				dreal mi = 0;
				for (integer y = 0;y < bins;++y)
				{
					dreal yMass = marginalHistogram(y, j);

					if (yMass > 0)
					{
						for (integer x = 0;x < bins;++x)
						{
							// We choose to do multiplications and division
							// instead of subtracting logarithms.
							// This way we hope to avoid possible cancellation
							// problems.

							dreal xMass = marginalHistogram(x, i);
							dreal xyMass = jointHistogram(x, y);

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

#endif
