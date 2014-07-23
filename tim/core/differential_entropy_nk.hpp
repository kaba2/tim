#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_HPP
#define TIM_DIFFERENTIAL_ENTROPY_NK_HPP

#include "tim/core/differential_entropy_nk.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/sequence_algorithms.h>
#include <pastel/sys/array_pointpolicy.h>

#include <pastel/math/matrix_tools.h>

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule.h>

#include <vector>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		typename SignalPtr_Range, 
		typename NormBijection>
	real differentialEntropyNk(
		const SignalPtr_Range& signalSet,
		const NormBijection& normBijection,
		integer* outIntrinsicDimension)
	{
		//const integer kNearest = 1;

		integer trials = signalSet.size();
		integer samples = minSamples(signalSet);
		integer dimension = signalSet.front()->dimension();

		const integer estimateSamples = samples * trials;

		// Generate codebook sizes.

		integer codebooks = dimension + 2;
		VectorD codebookSize(ofDimension(codebooks));
		for (integer i = 0;i < codebooks;++i)
		{
			// At least 10% of the samples must be in the
			// codebook, and at least 10% of the samples
			// must be outside the codebook.

			const real u = ((real)i / (codebooks - 1)) * 0.8 + 0.1;
			codebookSize[i] = (integer)(estimateSamples * u);
		}

		// Gather the point set.

		std::vector<const real*> pointSet;
		pointSet.reserve(estimateSamples);
		auto iter = signalSet.begin();
		auto iterEnd = signalSet.end();
		while(iter != iterEnd)
		{
			const Signal& signal = **iter;
			copy_n(
				signal.pointBegin(), samples,
				std::back_inserter(pointSet));
			
			++iter;
		}

		// Create a kd-tree.

		using Settings = PointKdTree_Settings<Pointer_Locator<real>>;
		typedef PointKdTree<Settings> KdTree;
		typedef KdTree::Point_ConstIterator Point_ConstIterator;
		typedef KdTree::Point Point;

		Pointer_Locator<real> locator(dimension);

		KdTree kdTree(locator);

		kdTree.insertRange(
			range(pointSet.begin(), pointSet.end()));
		kdTree.refine(SplitRule());

		// For each m, compute average log-distance alpha_m to the nearest 
		// codebook point for all points _not_ in the codebook 
		// (the nearest codebook point for a codebook point is the point 
		// itself).

		VectorD alphaSet(ofDimension(codebooks));
				
		for (integer m = 0;m < codebooks;++m)
		{
			// Select a random subset.

			integer subsetSize = codebookSize[m];
			randomSubset(
				pointSet.begin(), pointSet.end(),
				subsetSize);

			kdTree.erase();
			kdTree.insertRange(
				range(pointSet.begin(), pointSet.begin() + subsetSize));

			using Block = tbb::blocked_range<integer>;
			using Pair = std::pair<real, integer>;
			

			auto compute = [&](
				const Block& block,
				const Pair& start)
			{
				real alpha = start.first;
				integer acceptedSamples = start.second;
				for (integer i = block.begin(); i < block.end(); ++i)
				{
					real distance =
						searchNearest(
						kdTree,
						VectorD(ofDimension(dimension), withAliasing((real*)pointSet[i])),
						nullOutput(),
						allIndicator(),
						normBijection);

					// The logarithm of zero would give -infinity,
					// so we must avoid that. We remove all such
					// cases from the estimate.
					if (distance > 0)
					{
						alpha += normBijection.toLnNorm(distance);
						++acceptedSamples;
					}
				}

				return Pair(alpha, acceptedSamples);
			};

			auto reduce = [](const Pair& left, const Pair& right)
			{
				return Pair(
					left.first + right.first, 
					left.second + right.second);
			};

			// Find the distance to the nearest codebook point for
			// all points not in the codebook. Compute their
			// average logarithm.

			real alpha = 0;
			integer acceptedSamples = 0;

			std::tie(alpha, acceptedSamples) = 
				tbb::parallel_reduce(
					Block(subsetSize, estimateSamples),
					Pair(0, 0),
					compute,
					reduce);

			if (acceptedSamples > 0)
			{
				alpha /= acceptedSamples;
			}
			alphaSet[m] = alpha;
		}

		// Solve for the intrinsic dimensionality and
		// kappa.

		Matrix<real> e(2, codebooks);

		e.column(0) = -alphaSet;
		e.column(1) = 1;

		VectorD m = log(codebookSize);
		VectorD theta(ofDimension(2));

		real averageLogSize = sum(m) / codebooks;
		real averageAlpha = sum(alphaSet) / codebooks;

		// Find the integer dimensionality d that minimizes 
		// the cost function.

		real minCost = infinity<real>();
		integer d = 0;
		for (integer i = 0;i <= dimension;++i)
		{
			theta[0] = i;
			
			// Given a dimension i, this is the kappa
			// which minimizes the cost function.

			theta[1] = averageLogSize + i * averageAlpha;

			const real cost = dot(m - e * theta);
			if (cost < minCost)
			{
				d = i;
				minCost = cost;
			}
		}

		// Compute the differential entropy.

		real entropy = 0;

		if (d > 0)
		{
			entropy = 
				d * averageAlpha +
				averageLogSize + 
				constantEulerMascheroni<real>() +
				normBijection.lnVolumeUnitSphere(d);
		}
		else
		{
			entropy = -infinity<real>();
		}

		if (outIntrinsicDimension)
		{
			*outIntrinsicDimension = d;
		}

		return entropy;
	}

}

#endif
