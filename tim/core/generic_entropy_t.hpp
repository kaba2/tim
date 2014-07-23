#ifndef TIM_GENERIC_ENTROPY_T_HPP
#define TIM_GENERIC_ENTROPY_T_HPP

#include "tim/core/generic_entropy_t.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/counting_iterator.h>
#include <pastel/sys/range.h>

#include <pastel/geometry/pointkdtree.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <algorithm>
#include <numeric>

namespace Tim
{

	template <
		typename SignalPtr_Range, 
		typename EntropyAlgorithm,
		typename Real_Filter_Iterator>
	Signal temporalGenericEntropy(
		const SignalPtr_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE(odd(filter.size()));

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		// This function encapsulates the common
		// properties of the temporal entropy estimation 
		// algorithms based on k-nearest neighbors.

		if (signalSet.empty())
		{
			return Signal(0, 1);
		}

		Integer2 sharedTime = sharedTimeInterval(signalSet);
		integer estimateBegin = sharedTime[0];
		integer estimateEnd = sharedTime[1];
		integer samples = estimateEnd - estimateBegin;

		integer trials = signalSet.size();
		integer dimension = signalSet.front()->dimension();

		const integer totalSamples = samples * trials;

		ENSURE_OP(kNearest, <, totalSamples);

		// Copy the filter and replicate
		// the values to each trial.

		integer filterWidth = filter.size();
		integer filterRadius = filterWidth / 2;
		integer maxLocalFilterWidth = 
			std::min(filterWidth, samples);

		std::vector<real> copyFilter;

		copyFilter.reserve(filterWidth * trials);

		{
			Real_Filter_Iterator iter = filter.begin();
			Real_Filter_Iterator iterEnd = filter.end();
			while(iter != iterEnd)
			{
				std::fill_n(

					std::back_inserter(copyFilter), trials, *iter);
				++iter;
			}
		}

		Signal result(samples, 1, estimateBegin);
		integer missingValues = 0;

#pragma omp parallel
		{
		// Each worker thread has to create its own copy of
		// the signal point set. This is because the call
		// to SignalPointSet::setTimeWindow() is mutating.
		// This is a bit wasteful in memory, but I don't
		// know how else this could be done.

		Array<real> distanceArray(Vector2i(1, maxLocalFilterWidth * trials));

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
					distanceArray(i - searchBegin) =
						searchNearest(
						pointSet.kdTree(),
						*(pointSet.begin() + i),
						nullOutput(),
						predicateIndicator(*(pointSet.begin() + i), NotEqualTo()),
						entropyAlgorithm.normBijection()).
						kNearest(kNearest);
				}
			};

			tbb::parallel_for(
				Block(searchBegin, searchEnd),
				search);

			// After we have found the distances, we simply evaluate
			// the generic entropy estimator over the samples of
			// the current time instant.

			real weightSum = 0;
			const integer filterOffset = tFilterOffset * trials;
			real estimate = 0;
			for (integer i = 0;i < windowSamples;++i)
			{
				// Points that are at identical positions do not
				// provide any information. Such samples are
				// not taken in the estimate.
				if (distanceArray(i) > 0)
				{
					real weight = copyFilter[i + filterOffset];

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

				result.data()(t - estimateBegin) = nan<real>();
				++missingValues;
			}
		}
		}

		// Reconstruct the NaN's.

		reconstruct(
			range(result.data().begin(), result.data().end()));

		return result;
	}

	template <
		typename SignalPtr_Range, 
		typename EntropyAlgorithm>
	Signal temporalGenericEntropy(
		const SignalPtr_Range& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		integer kNearest)
	{
		return Tim::temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
