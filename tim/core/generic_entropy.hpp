#ifndef TIM_GENERIC_ENTROPY_HPP
#define TIM_GENERIC_ENTROPY_HPP

#include "tim/core/generic_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/iterator/constant_iterator.h>
#include <pastel/sys/iterator/counting_iterator.h>
#include <pastel/sys/range.h>
#include <pastel/sys/indicator/predicate_indicator.h>

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest_kdtree.h>

#include <algorithm>
#include <numeric>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		typename Signal_Range,
		typename EntropyAlgorithm>
	real genericEntropy(
		const Signal_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (signalSet.empty())
		{
			return (real)Nan();
		}

		// This function encapsulates the common
		// properties of the entropy estimation 
		// algorithms based on k-nearest neighbors.

		SignalPointSet pointSet(signalSet);

		integer trials = signalSet.size();
		integer samples = pointSet.samples();
		integer dimension = signalSet.front()->dimension();

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
		using Pair = std::pair<real, integer>;
		
		auto compute = [&](
			const Block& block,
			const Pair& start)
		{
			real estimate = start.first;
			integer acceptedSamples = start.second;

			for (integer i = block.begin();i < block.end();++i)
			{
				auto query = indexedPointSet[i];

				Vector<real> queryPoint(
					ofDimension(pointSet.dimension()),
					withAliasing((real*)(query->point())));

				// Find the distance to the k:th nearest neighbor.
				real distance2 =
					searchNearest(
						pointSet.kdTree(),
						queryPoint,
						PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
						PASTEL_TAG(normBijection), entropyAlgorithm.normBijection(),
						PASTEL_TAG(kNearest), kNearest
					).first;

				// Points that are at identical positions do not
				// provide any information. Such samples are
				// not taken in the estimate.
				if (distance2 > 0)
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

		real estimate = 0;
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
			estimate = (real)Nan();
		}

		return estimate;
	}

}

#endif
