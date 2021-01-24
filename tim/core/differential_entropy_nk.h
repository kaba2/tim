// Description: Differential entropy estimation
// Detail: Nilsson-Kleijn manifold nearest neighbor estimator

#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_H
#define TIM_DIFFERENTIAL_ENTROPY_NK_H

#include "tim/core/mytypes.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/range.h>
#include <pastel/sys/sequence/sequence_algorithms.h>
#include <pastel/sys/locator/pointer_locator.h>

#include <pastel/math/matrix/matrix_tools.h>

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/splitrule/slidingmidpoint_splitrule.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <vector>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		ranges::forward_range Signal_Range, 
		typename Norm>
	dreal differentialEntropyNk(
		const Signal_Range& signalSet,
		const Norm& norm,
		integer* outIntrinsicDimension = 0)
	{
		//const integer kNearest = 1;

		integer trials = ranges::size(signalSet);
		integer samples = minSamples(signalSet);
		integer dimension = std::begin(signalSet)->dimension();

		const integer estimateSamples = samples * trials;

		// Generate codebook sizes.

		integer codebooks = dimension + 2;
		VectorD codebookSize(ofDimension(codebooks));
		for (integer i = 0;i < codebooks;++i)
		{
			// At least 10% of the samples must be in the
			// codebook, and at least 10% of the samples
			// must be outside the codebook.

			const dreal u = ((dreal)i / (codebooks - 1)) * 0.8 + 0.1;
			codebookSize[i] = (integer)(estimateSamples * u);
		}

		// Gather the point set.

		std::vector<const dreal*> pointSet;
		pointSet.reserve(estimateSamples);
		auto iter = ranges::begin(signalSet);
		auto iterEnd = signalSet.end();
		while(iter != iterEnd)
		{
			const Signal& signal = *iter;
			copy_n(
				std::begin(signal.pointRange()), samples,
				std::back_inserter(pointSet));
			
			++iter;
		}

		// Create a kd-tree.

		using Settings = PointKdTree_Settings<Pointer_Locator<dreal>>;
		typedef PointKdTree<Settings> KdTree;
		typedef KdTree::Point_ConstIterator Point_ConstIterator;
		typedef KdTree::Point Point;

		Pointer_Locator<dreal> locator(dimension);

		KdTree kdTree(locator);

		kdTree.insertSet(pointSet);
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
			kdTree.insertSet(
				range(pointSet.begin(), pointSet.begin() + subsetSize));

			using Block = tbb::blocked_range<integer>;
			using Pair = std::pair<dreal, integer>;
			

			auto compute = [&](
				const Block& block,
				const Pair& start)
			{
				dreal alpha = start.first;
				integer acceptedSamples = start.second;
				for (integer i = block.begin(); i < block.end(); ++i)
				{
					VectorD queryPoint(
						ofDimension(dimension),
						withAliasing((dreal*)pointSet[i]));

					auto distance =
						searchNearest(
							kdTreeNearestSet(kdTree),
							queryPoint,
							PASTEL_TAG(norm), norm
						).first;

					// The logarithm of zero would give -infinity,
					// so we must avoid that. We remove all such
					// cases from the estimate.
					if ((dreal)distance > 0)
					{
						alpha += std::log((dreal)distance);
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

			dreal alpha = 0;
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

		Matrix<dreal> e(codebooks, 2);

		e.array().col(0) = -asColumnMatrix(alphaSet);
		e.array().col(1) = 1;

		VectorD m = log(codebookSize);
		VectorD theta(ofDimension(2));

		dreal averageLogSize = sum(m) / codebooks;
		dreal averageAlpha = sum(alphaSet) / codebooks;

		// Find the integer dimensionality d that minimizes 
		// the cost function.

		dreal minCost = Infinity();
		integer d = 0;
		for (integer i = 0;i <= dimension;++i)
		{
			theta[0] = i;
			
			// Given a dimension i, this is the kappa
			// which minimizes the cost function.

			theta[1] = averageLogSize + i * averageAlpha;

			const dreal cost = dot(m - e * theta);
			if (cost < minCost)
			{
				d = i;
				minCost = cost;
			}
		}

		// Compute the differential entropy.

		dreal entropy = 0;

		if (d > 0)
		{
			entropy = 
				d * averageAlpha +
				averageLogSize + 
				constantEulerMascheroni<dreal>() +
				lnVolumeUnitSphere(norm, d);
		}
		else
		{
			entropy = -Infinity();
		}

		if (outIntrinsicDimension)
		{
			*outIntrinsicDimension = d;
		}

		return entropy;
	}

}

#endif
