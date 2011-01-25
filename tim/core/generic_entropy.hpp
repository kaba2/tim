#ifndef TIM_GENERIC_ENTROPY_HPP
#define TIM_GENERIC_ENTROPY_HPP

#include "tim/core/generic_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/countingiterator.h>
#include <pastel/sys/iterator_range.h>

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

#include <algorithm>
#include <numeric>

namespace Tim
{

	template <
		typename SignalPtr_Iterator,
		typename EntropyAlgorithm>
	real genericEntropy(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet,
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

		// This is done to avoid parallelization
		// issues with iterator range caching.

		signalSet.updateCache();

		SignalPointSet pointSet(signalSet, true);

		const integer trials = signalSet.size();
		const integer samples = pointSet.samples();
		const integer dimension = signalSet.front()->dimension();
		const integer estimateSamples = samples * trials;

		Array<real> distanceArray(1, estimateSamples);

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

		integer acceptedSamples = 0;
		real estimate = 0;
#pragma omp parallel for reduction(+ : estimate, acceptedSamples)
		for (integer i = 0;i < estimateSamples;++i)
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
