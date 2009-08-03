#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_HPP
#define TIM_MUTUAL_INFORMATION_KRASKOW_HPP

#include "tim/core/mutual_information_kraskow.h"

#include <pastel/sys/constantiterator.h>

namespace Tim
{

	template <
		typename Signal_A_ForwardIterator,
		typename Signal_B_ForwardIterator,
		typename Lag_ForwardIterator>
	real mutualInformation(
		const Signal_A_ForwardIterator& aTrialBegin,
		const Signal_B_ForwardIterator& bTrialBegin,
		integer trials,
		const Lag_ForwardIterator& lagBegin,
		integer lags,
		integer sigma,
		integer kNearest)
	{
		ENSURE_OP(trials, >, 0);
		ENSURE_OP(lags, >, 0);
		ENSURE_OP(kNearest, >, 0);
		PENSURE(nEqualDimension(aTrialBegin, trials));
		PENSURE(nEqualDimension(bTrialBegin, trials));

		if (trials == 0)
		{
			return 0;
		}

		const integer samples = std::min(
			minSamples(aTrialBegin, trials),
			minSamples(bTrialBegin, trials));

		if (sigma < 0)
		{
			sigma = samples;
		}

		const integer maxPointsPerNode = 16;
		SlidingMidpoint2_SplitRule splitRule;

		const Infinity_NormBijection<real> normBijection;
		Array<2, real> distanceArray(1, samples);

		typedef PointKdTree<N, Real, 
			Detail_AllNearestNeighborsKdTree::PointListPolicy<N, Real> > KdTree;
		typedef typename KdTree::ConstObjectIterator ConstTreeIterator;
		typedef typename KdTree::Object Object;
		typedef CountingIterator<const Point<N, Real>*> SequenceIterator;

		Signal_A_ForwardIterator aIter = aTrialBegin;
		Signal_B_ForwardIterator bIter = bTrialBegin;
		Lag_Iterator lagIter = lagBegin;
		for (integer j = 0;j < lags;++j)
		{
			std::vector<SignalPtr> jointSet;
			jointSet.reserve(trials);

			integer bLag = *lagIter;
			for (integer i = 0;i < trials;++i)
			{
				jointSet.push_back(merge(*aIter, *bIter, bLag));
				++aIter;
				++bIter;
			}

			// Create the kd-tree subdivision based on all samples.

			std::vector<PointD> pointSet;
			for (integer i = 0;i < trials;++i)
			{
				constructPointSet(jointSet[i], pointSet);
			}

			KdTree kdTree(maxPointsPerNode, dimension);

			kdTree.insert(
				SequenceIterator(&pointSet[0]), 
				SequenceIterator(&pointSet[0] + pointSet.size()));

			kdTree.refine(128, maxPointsPerNode, splitRule);

			searchAllNeighborsKdTree(
				kdTree,
				DepthFirst_SearchAlgorithm_PointKdTree(),
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				kNearest - 1,
				kNearest,
				infinity<real>(),
				maxRelativeError,
				normBijection,
				0,
				&distanceArray);

			real estimate = 0;

			std::vector<integer> countSet(samples, 0);
			integer dimensionOffset = 0;
			integer signals = 0;

			Signal_ForwardIterator signalIter = signalBegin;
			while(signalIter != signalEnd)
			{
				const integer dimension = (*signalIter)->dimension();
				
				constructPointSet(
					jointSignal,
					0, samples,
					dimensionOffset, 
					dimensionOffset + dimension,
					pointSet);

				countAllNeighborsKdTree(
					pointSet,
					CountingIterator<integer>(0),
					CountingIterator<integer>(samples),
					distanceArray.begin(),
					normBijection,
					16,
					countSet.begin());

	#pragma omp parallel for reduction(+ : estimate)
				for (integer j = 0;j < samples;++j)
				{
					estimate -= digamma<real>(countSet[j]);
				}

				dimensionOffset += dimension;
				++signalIter;
				++signals;
			}

			estimate /= samples;
			estimate += digamma<real>(kNearest);
			estimate += (signals - 1) * digamma<real>(samples);

			++lagIter;
		}

		return estimate;
	}

}

#endif
