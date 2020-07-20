// Description: Temporal generic entropy estimation
// Detail: Encapsulates properties common to k-nn entropy estimators.

#ifndef TIM_GENERIC_ENTROPY_T_H
#define TIM_GENERIC_ENTROPY_T_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <numeric>

namespace Tim
{

	//! Computes temporal generic entropy of a signal.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0

	signalSet:
	An ensemble of signals representing trials
	of the same experiment.

	entropyAlgorithm:
	Encapsulates the specifics of the used
	entropy algorithm.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate generic entropy.

	filter:
	An array of coefficients by which to weight the results
	in the time-window. The center of the array corresponds 
	to the current time instant. The width of the array can 
	be arbitrary but must be odd. The coefficients must sum
	to a non-zero value.
	*/
	template <
		ranges::forward_range Signal_Range, 
		typename EntropyAlgorithm,
		ranges::forward_range Filter_Range>
	SignalData temporalGenericEntropy(
		const Signal_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		integer kNearest,
		const Filter_Range& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE(odd(ranges::size(filter)));

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		// This function encapsulates the common
		// properties of the temporal entropy estimation 
		// algorithms based on k-nearest neighbors.

		if (ranges::empty(signalSet))
		{
			return SignalData(0, 1);
		}

		auto norm = entropyAlgorithm.norm();
		using Distance = decltype(norm());

		Integer2 sharedTime = sharedTimeInterval(signalSet);
		integer estimateBegin = sharedTime[0];
		integer estimateEnd = sharedTime[1];
		integer samples = estimateEnd - estimateBegin;

		integer trials = ranges::size(signalSet);
		integer dimension = std::begin(signalSet)->dimension();

		const integer totalSamples = samples * trials;

		ENSURE_OP(kNearest, <, totalSamples);

		// Copy the filter and replicate
		// the values to each trial.

		integer filterWidth = ranges::size(filter);
		integer filterRadius = filterWidth / 2;
		integer maxLocalFilterWidth = 
			std::min(filterWidth, samples);

		std::vector<dreal> copyFilter;

		copyFilter.reserve(filterWidth * trials);

		{
			auto iter = std::begin(filter);
			auto iterEnd = std::end(filter);
			while(iter != iterEnd)
			{
				std::fill_n(
					std::back_inserter(copyFilter), trials, *iter);
				++iter;
			}
		}

		SignalData result(samples, 1, estimateBegin);
		integer missingValues = 0;

#pragma omp parallel
		{
		// Each worker thread has to create its own copy of
		// the signal point set. This is because the call
		// to SignalPointSet::setTimeWindow() is mutating.
		// This is a bit wasteful in memory, but I don't
		// know how else this could be done.

		Array<Distance> distanceArray(Vector2i(1, maxLocalFilterWidth * trials));

		SignalPointSet pointSet(signalSet);

#pragma omp for reduction(+ : missingValues)
		for (integer t = estimateBegin;t < estimateEnd;++t)
		{
			// Update the position of the time-window.

			pointSet.setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);

			integer tBegin = pointSet.windowBegin();
			integer tEnd = pointSet.windowEnd();
			integer tWidth = tEnd - tBegin;
			integer tLocalFilterBegin = std::max(t - filterRadius, tBegin) - tBegin;
			integer tLocalFilterEnd = std::min(t + filterRadius + 1, tEnd) - tBegin;
			integer tFilterDelta = tBegin - (t - filterRadius);
			integer tFilterOffset = std::max(tFilterDelta, (integer)0);

			const integer windowSamples = (tLocalFilterEnd - tLocalFilterBegin) * trials;
			
			// For each point at the current time instant in all
			// ensemble signals, find the distance to the k:th nearest 
			// neighbor. Note that the SignalPointSet stores the
			// point iterators interleaved so that for a given time instant
			// the samples of ensemble signals are listed sequentially.
			// I.e. if the ensemble signals are A, B and C, then
			// SignalPointSet stores point iterators to 
			// A(1), B(1), C(1), A(2), B(2), C(2), etc.
			// That is, the distance between subsequent samples of a
			// specific signal are 'trials' samples away.

			using Block = tbb::blocked_range<integer>;

			integer searchBegin = tLocalFilterBegin * trials;
			integer searchEnd = tLocalFilterEnd * trials;

			auto search = [&](const Block& block)
			{
				for (integer i = block.begin(); i < block.end(); ++i)
				{
					auto query = *(pointSet.begin() + i);

					Vector<dreal> queryPoint(
						ofDimension(pointSet.dimension()),
						withAliasing((dreal*)(query->point())));

					distanceArray(i - searchBegin) =
						searchNearest(
							kdTreeNearestSet(pointSet.kdTree()),
							queryPoint,
							PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
							PASTEL_TAG(norm), entropyAlgorithm.norm(),
							PASTEL_TAG(kNearest), kNearest
						).first;
				}
			};

			tbb::parallel_for(
				Block(searchBegin, searchEnd),
				search);

			// After we have found the distances, we simply evaluate
			// the generic entropy estimator over the samples of
			// the current time instant.

			dreal weightSum = 0;
			const integer filterOffset = tFilterOffset * trials;
			dreal estimate = 0;
			for (integer i = 0;i < windowSamples;++i)
			{
				// Points that are at identical positions do not
				// provide any information. Such samples are
				// not taken in the estimate.
				if ((dreal)distanceArray(i) > 0)
				{
					dreal weight = copyFilter[i + filterOffset];

					estimate += weight * entropyAlgorithm.sumTerm(distanceArray(i));
					weightSum += weight;
				}
			}
			if (weightSum != 0)
			{
				result.data()(t - estimateBegin) = 
					entropyAlgorithm.finishEstimate(
					estimate / weightSum, dimension, 
					kNearest, tWidth * trials);
			}
			else
			{
				// If all distances were zero, we can't say
				// anything about generic entropy. This is
				// marked with a NaN. We will later attempt
				// to reconstruct these values.

				result.data()(t - estimateBegin) = (dreal)Nan();
				++missingValues;
			}
		}
		}

		// Reconstruct the NaN's.

		reconstruct(result.data().range());

		return result;
	}

	//! Computes temporal generic entropy of a signal.
	/*!
	This is a convenience function that calls:

	temporalGenericEntropy(
		signalSet,
		entropyAlgorithm,
		timeWindowRadius,
		result,
		kNearest,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
	template <
		ranges::forward_range Signal_Range, 
		typename EntropyAlgorithm>
	SignalData temporalGenericEntropy(
		const Signal_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		integer kNearest = 1)
	{
		return Tim::temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			constantRange((dreal)1, 1));
	}

}

#endif
