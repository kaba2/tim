#ifndef TIM_ENTROPY_COMBINATION2_HPP
#define TIM_ENTROPY_COMBINATION2_HPP

#include "tim/core/entropy_combination2.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/forwardrange.h>

#include <numeric>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	real entropyCombination2(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(lagSet.size(), ==, signalSet.height());

		if (signalSet.empty() || rangeSet.empty())
		{
			return 0;
		}

		// Construct the joint signal.

		const integer trials = signalSet.width();

		std::vector<SignalPtr> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, std::back_inserter(jointSignalSet), lagSet);

		const integer samples = jointSignalSet.front()->samples();
		if (samples == 0)
		{
			return 0;
		}

		const integer signals = signalSet.height();

		// Find out the dimension ranges of the marginal
		// signals.

		std::vector<integer> offsetSet;
		offsetSet.reserve(signals + 1);
		offsetSet.push_back(0);
		for (integer i = 1;i < signals + 1;++i)
		{
			offsetSet.push_back(offsetSet[i - 1] + signalSet(0, i - 1)->dimension());
		}

		const integer estimateSamples = samples * trials;
		const integer marginals = rangeSet.size();

		// Construct point sets

		SignalPointSet jointPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end()), true);
	
		std::vector<Integer3> copyRangeSet(
			rangeSet.begin(), rangeSet.end());

		std::vector<SignalPointSetPtr> pointSet;
		pointSet.reserve(marginals);
		std::vector<integer> macroDimensionSet;
		macroDimensionSet.reserve(marginals);

		real sumWeight = 0;
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];
			pointSet.push_back(
				SignalPointSetPtr(new SignalPointSet(
				forwardRange(jointSignalSet.begin(), jointSignalSet.end()), 
				true,
				offsetSet[range[0]], offsetSet[range[1]])));
			macroDimensionSet.push_back(range[1] - range[0]);

			sumWeight += range[2];
		}

		// It is essential that the used norm is the
		// infinity norm.

		Infinity_NormBijection<real> normBijection;

		// Start estimation.

		typedef SignalPointSet::ConstObjectIterator ConstObjectIterator;

		Array<ConstObjectIterator, 2> nearestArray(1, estimateSamples);

		searchAllNeighbors(
			jointPointSet.kdTree(),
			DepthFirst_SearchAlgorithm_PointKdTree(),
			randomAccessRange(jointPointSet.begin(), jointPointSet.end()),
			kNearest - 1,
			kNearest, 
			randomAccessRange(constantIterator(infinity<real>()), estimateSamples),
			0,
			normBijection,
			&nearestArray);

		std::vector<real> distanceArray(estimateSamples);

		real estimate = 0;
		std::vector<integer> countSet(estimateSamples, 0);
		for (integer i = 0;i < marginals;++i)
		{
			const integer marginalDimension = pointSet[i]->dimension();
			const integer marginalOffset = pointSet[i]->dimensionBegin();

#pragma omp parallel for
			for (integer j = 0;j < estimateSamples;++j)
			{
				const ConstObjectIterator from = jointPointSet.begin()[j];
				const ConstObjectIterator to = nearestArray(j);
				distanceArray[j] = distance2(
					from->object() + marginalOffset,
					to->object() + marginalOffset,
					marginalDimension,
					normBijection);
			}

			countAllNeighbors(
				pointSet[i]->kdTree(),
				randomAccessRange(pointSet[i]->begin(), pointSet[i]->end()),
				randomAccessRange(distanceArray.begin(), estimateSamples),
				normBijection,
				countSet.begin());

			const integer macroDimension = macroDimensionSet[i];

			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate)
			for (integer j = 0;j < estimateSamples;++j)
			{
				const integer k = std::max(countSet[j] - 1, 1);

				signalEstimate += digamma<real>(k);
				//signalEstimate -= (real)(macroDimension - 1) / k;
			}

			estimate -= signalEstimate * copyRangeSet[i][2];
		}
		estimate /= estimateSamples;
		estimate += digamma<real>(kNearest);
		//estimate -= (real)(marginals - 1) / kNearest;
		estimate -= 1;

		estimate += (sumWeight - 1) * digamma<real>(estimateSamples);

		return estimate;
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	real entropyCombination2(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet)
	{
		return Tim::entropyCombination2(
			signalSet,
			rangeSet,
			forwardRange(constantIterator(0), signalSet.height()),
			kNearest);
	}

}

#endif
