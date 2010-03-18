#ifndef TIM_GENERIC_ENTROPY_T_HPP
#define TIM_GENERIC_ENTROPY_T_HPP

#include "tim/core/generic_entropy_t.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/countingiterator.h>
#include <pastel/sys/randomaccessrange.h>

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

#include <algorithm>
#include <numeric>

namespace Tim
{

	template <
		typename SignalPtr_Iterator, 
		typename EntropyAlgorithm,
		typename Real_OutputIterator,
		typename Real_Filter_Iterator>
	integer temporalGenericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE(odd(filter.size()));

		// This function encapsulates the common
		// properties of the temporal entropy estimation 
		// algorithms based on k-nearest neighbors.

		if (signalSet.empty())
		{
			return 0;
		}

		// This is done to avoid parallelization
		// issues with iterator range caching.

		signalSet.updateCache();

		const integer trials = signalSet.size();
		const integer samples = minSamples(signalSet);
		const integer dimension = signalSet.front()->dimension();
		const integer totalSamples = samples * trials;

		ENSURE_OP(kNearest, <, totalSamples);

		// Copy the filter and replicate
		// the values to each trial.

		const integer filterWidth = filter.size();
		const integer filterRadius = filterWidth / 2;
		const integer maxLocalFilterWidth = 
			std::min(filterWidth, samples);

		std::vector<real> copyFilter;
		copyFilter.reserve(filterWidth * trials);

		{
			Real_Filter_Iterator iter = filter.begin();
			const Real_Filter_Iterator iterEnd = filter.end();
			while(iter != iterEnd)
			{
				std::fill_n(
					std::back_inserter(copyFilter), trials, *iter);
				++iter;
			}
		}

		// We create an own array to hold the results since the
		// 'result' iterator is not necessarily random-access
		// (which is needed for parallelization below).

		std::vector<real> estimateSet(samples);
		integer missingValues = 0;

#pragma omp parallel
		{
		// Each worker thread has to create its own copy of
		// the signal point set. This is because the call
		// to SignalPointSet::setTimeWindow() is mutating.
		// This is a bit wasteful in memory, but I don't
		// know how else this could be done.

		Array<real, 2> distanceArray(1, maxLocalFilterWidth * trials);

		SignalPointSet pointSet(signalSet);

#pragma omp for reduction(+ : missingValues)
		for (integer t = 0;t < samples;++t)
		{
			// Update the position of the time-window.

			pointSet.setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);

			const integer tBegin = pointSet.timeBegin();
			const integer tEnd = pointSet.timeEnd();
			const integer tWidth = tEnd - tBegin;
			const integer tLocalFilterBegin = std::max(t - filterRadius, tBegin) - tBegin;
			const integer tLocalFilterEnd = std::min(t + filterRadius + 1, tEnd) - tBegin;
			const integer tFilterDelta = tBegin - (t - filterRadius);
			const integer tFilterOffset = std::max(tFilterDelta, 0);
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

			searchAllNeighbors(
				pointSet.kdTree(),
				randomAccessRange(
				pointSet.begin() + tLocalFilterBegin * trials, 
				pointSet.begin() + tLocalFilterEnd * trials),
				kNearest - 1,
				kNearest, 
				0,
				&distanceArray,
				constantRange(infinity<real>(), windowSamples),
				0,
				entropyAlgorithm.normBijection());

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
					const real weight = copyFilter[i + filterOffset];
					estimate += weight * entropyAlgorithm.sumTerm(distanceArray(i));
					weightSum += weight;
				}
			}
			if (weightSum != 0)
			{
				estimateSet[t] = entropyAlgorithm.finishEstimate(
					estimate / weightSum, dimension, 
					kNearest, tWidth * trials);
			}
			else
			{
				// If all distances were zero, we can't say
				// anything about generic entropy. This is
				// marked with a NaN. We will later attempt
				// to reconstruct these values.

				estimateSet[t] = nan<real>();
				++missingValues;
			}
		}
		}

		// Reconstruct the NaN's.

		reconstruct(
			forwardRange(estimateSet.begin(), estimateSet.end()));

		// Copy the results to the output.

		std::copy(estimateSet.begin(), estimateSet.end(), result);

		return missingValues;
	}

	template <
		typename SignalPtr_Iterator, 
		typename EntropyAlgorithm,
		typename Real_OutputIterator>
	integer temporalGenericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		return Tim::temporalGenericEntropy(
			signalSet,
			entropyAlgorithm,
			timeWindowRadius,
			result,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
