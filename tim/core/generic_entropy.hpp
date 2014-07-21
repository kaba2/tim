#ifndef TIM_GENERIC_ENTROPY_HPP
#define TIM_GENERIC_ENTROPY_HPP

#include "tim/core/generic_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/counting_iterator.h>
#include <pastel/sys/range.h>
#include <pastel/sys/predicate_indicator.h>

#include <pastel/geometry/pointkdtree.h>

#include <algorithm>
#include <numeric>

#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>

namespace Tim
{

	template <
		typename SignalPtr_Iterator,
		typename EntropyAlgorithm>
	real genericEntropy(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (signalSet.empty())
		{
			return nan<real>();
		}

		// This function encapsulates the common
		// properties of the entropy estimation 
		// algorithms based on k-nearest neighbors.

		SignalPointSet pointSet(signalSet);

		const integer trials = signalSet.size();
		const integer samples = pointSet.samples();
		const integer dimension = signalSet.front()->dimension();
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
				// Find the distance to the k:th nearest neighbor.
				real distance2 =
					searchNearest(
					pointSet.kdTree(),
					indexedPointSet[i],
					nullOutput(),
					predicateIndicator(indexedPointSet[i], NotEqualTo()),
					entropyAlgorithm.normBijection())
					.kNearest(kNearest);

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
			estimate = nan<real>();
		}

		return estimate;
	}

}

#endif
