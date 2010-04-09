#ifndef TIM_DIFFERENTIAL_ENTROPY_NK_HPP
#define TIM_DIFFERENTIAL_ENTROPY_NK_HPP

#include "tim/core/differential_entropy_nk.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/stdext_algorithm.h>

#include <pastel/math/matrix_tools.h>

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>
#include <pastel/geometry/array_pointpolicy.h>

#include <vector>

namespace Tim
{

	template <
		typename SignalPtr_Iterator, 
		typename NormBijection>
	real differentialEntropyNk(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const NormBijection& normBijection,
		integer* outIntrinsicDimension)
	{
		const integer kNearest = 1;

		const integer trials = signalSet.size();
		const integer samples = minSamples(signalSet);
		const integer dimension = signalSet.front()->dimension();
		const integer estimateSamples = samples * trials;

		// Generate codebook sizes.

		const integer codebooks = dimension + 2;
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
		SignalPtr_Iterator iter = signalSet.begin();
		const SignalPtr_Iterator iterEnd = signalSet.end();
		while(iter != iterEnd)
		{
			const SignalPtr& signal = *iter;
			StdExt::copy_n(
				signal->pointBegin(), samples,
				std::back_inserter(pointSet));
			
			++iter;
		}

		// Create a kd-tree.

		typedef PointKdTree<real, Dynamic, Array_PointPolicy<real> > KdTree;
		typedef KdTree::ConstObjectIterator ConstObjectIterator;
		typedef KdTree::Object Object;

		KdTree kdTree(ofDimension(dimension));

		kdTree.insert(pointSet.begin(), pointSet.end());
		kdTree.refine(SplitRule());

		// For each m, compute average log-distance alpha_m to the nearest 
		// codebook point for all points _not_ in the codebook 
		// (the nearest codebook point for a codebook point is the point 
		// itself).

		VectorD alphaSet(ofDimension(codebooks));
				
		for (integer m = 0;m < codebooks;++m)
		{
			// Select a random subset.

			const integer subsetSize = codebookSize[m];
			StdExt::subset(
				pointSet.begin(), pointSet.end(),
				subsetSize);

			kdTree.eraseObjects();
			kdTree.insert(pointSet.begin(), pointSet.begin() + subsetSize);
			
			// Find the distance to the nearest codebook point for
			// all points not in the codebook. Compute their
			// average logarithm.

			integer acceptedSamples = 0;
			real alpha = 0;
#pragma omp parallel for reduction(+ : alpha, acceptedSamples)
			for (integer i = subsetSize;i < estimateSamples;++i)
			{
				const real distance =
					searchNearestOne(
					kdTree,
					VectorD(ofDimension(dimension), withAliasing((real*)pointSet[i])),
					infinity<real>(), 0,
					Always_AcceptPoint<ConstObjectIterator>(),
					8,
					normBijection).key();

				// The logarithm of zero would give -infinity,
				// so we must avoid that. We remove all such
				// cases from the estimate.
				if (distance > 0)
				{
					alpha += normBijection.toLnNorm(distance);
					++acceptedSamples;
				}
			}
			if (acceptedSamples > 0)
			{
				alpha /= acceptedSamples;
			}
			alphaSet[m] = alpha;
		}

		// Solve for the intrinsic dimensionality and
		// kappa.

		MatrixD e(2, codebooks);

		e[0] = -alphaSet;
		e[1] = 1;

		const VectorD m = log(codebookSize);
		VectorD theta(ofDimension(2));

		const real averageLogSize = sum(m) / codebooks;
		const real averageAlpha = sum(alphaSet) / codebooks;

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

			const real cost = dot(m - theta * e);
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
