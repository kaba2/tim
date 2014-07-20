#ifndef TIM_GENERIC_ENTROPY_HPP
#define TIM_GENERIC_ENTROPY_HPP

#include "tim/core/generic_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/counting_iterator.h>
#include <pastel/sys/range.h>

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

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

		Array<real> distanceArray(Vector2i(1, estimateSamples));

		// Find the distance to the k:th nearest neighbor for all points.

		searchAllNeighbors(
			pointSet.kdTree(),
			range(pointSet.begin(), pointSet.end()),
			kNearest - 1,
			kNearest, 
			(Array<Point_ConstIterator>*)0,
			&distanceArray,
			constantRange(infinity<real>(), estimateSamples),
			0,
			entropyAlgorithm.normBijection());

		// After we have found the distances, we simply evaluate
		// the generic entropy estimator over all samples.

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
				// Points that are at identical positions do not
				// provide any information. Such samples are
				// not taken in the estimate.
				if (distanceArray(i) > 0)
				{
					estimate += entropyAlgorithm.sumTerm(distanceArray(i));
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
