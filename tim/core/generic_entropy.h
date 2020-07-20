// Description: Generic entropy estimation
// Detail: Encapsulates properties common to k-nn entropy estimators.

#ifndef TIM_GENERIC_ENTROPY_H
#define TIM_GENERIC_ENTROPY_H

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/range.h>
#include <pastel/sys/indicator/predicate_indicator.h>

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <algorithm>
#include <numeric>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	//! Generic entropy of a signal.
	/*!
	Preconditions:
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	entropyAlgorithm:
	Encapsulates the specifics of the used
	entropy algorithm.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate generic entropy.

	Returns:
	A generic entropy estimate if successful,
	NaN otherwise. The estimation may fail only
	if all points are at the same position or
	there are no samples to estimate from.
	*/
	template <
		ranges::forward_range Signal_Range,
		typename EntropyAlgorithm>
	dreal genericEntropy(
		const Signal_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer kNearest = 1)
	{
		ENSURE_OP(kNearest, >, 0);

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (ranges::empty(signalSet))
		{
			return (dreal)Nan();
		}

		auto norm = entropyAlgorithm.norm();
		using Distance = decltype(norm());

		// This function encapsulates the common
		// properties of the entropy estimation 
		// algorithms based on k-nearest neighbors.

		SignalPointSet pointSet(signalSet);

		integer trials = ranges::size(signalSet);
		integer samples = pointSet.samples();
		integer dimension = std::begin(signalSet)->dimension();

		const integer estimateSamples = samples * trials;

		// Store the point-iterators into an array
		// for random-access for parallel_for.
		std::vector<Point_ConstIterator> indexedPointSet;
		indexedPointSet.reserve(pointSet.kdTree().points());
		for (auto i = pointSet.kdTree().begin(); i != pointSet.kdTree().end(); ++i)
		{
			indexedPointSet.emplace_back(i);
		}

		using Block = tbb::blocked_range<integer>;
		using Pair = std::pair<dreal, integer>;
		
		auto compute = [&](
			const Block& block,
			const Pair& start)
		{
			dreal estimate = start.first;
			integer acceptedSamples = start.second;

			for (integer i = block.begin();i < block.end();++i)
			{
				auto query = indexedPointSet[i];

				Vector<dreal> queryPoint(
					ofDimension(pointSet.dimension()),
					withAliasing((dreal*)(query->point())));

				// Find the distance to the k:th nearest neighbor.
				Distance distance2 =
					searchNearest(
						kdTreeNearestSet(pointSet.kdTree()),
						queryPoint,
						PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
						PASTEL_TAG(norm), entropyAlgorithm.norm(),
						PASTEL_TAG(kNearest), kNearest
					).first;

				// Points that are at identical positions do not
				// provide any information. Such samples are
				// not taken in the estimate.
				if ((dreal)distance2 > 0)
				{
					estimate += entropyAlgorithm.sumTerm(distance2);
					++acceptedSamples;
				}
			}

			return Pair(estimate, acceptedSamples);
		};

		auto reduce = [](const Pair& left, const Pair& right)
		{
			return Pair(
				left.first + right.first, 
				left.second + right.second);
		};

		dreal estimate = 0;
		integer acceptedSamples = 0;

		std::tie(estimate, acceptedSamples) = 
			tbb::parallel_reduce(
				Block(0, estimateSamples),
				Pair(0, 0),
				compute,
				reduce);

		if (acceptedSamples > 0)
		{
			estimate = entropyAlgorithm.finishEstimate(
				estimate / acceptedSamples, dimension, kNearest, 
				estimateSamples);
		}
		else
		{
			// If all distances were zero, we can't say
			// anything about generic entropy. This is
			// marked with a NaN.
			estimate = (dreal)Nan();
		}

		return estimate;
	}

}

#endif
