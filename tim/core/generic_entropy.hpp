#ifndef TIM_GENERIC_ENTROPY_HPP
#define TIM_GENERIC_ENTROPY_HPP

#include "tim/core/generic_entropy.h"
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

	// Temporal generic entropy
	// ------------------------

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
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);

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

		Array<real, 2> distanceArray(1, trials);
		SignalPointSet pointSet(signalSet);

#pragma omp for reduction(+ : missingValues)
		for (integer t = 0;t < samples;++t)
		{
			// Update the position of the time-window.

			pointSet.setTimeWindow(t - timeWindowRadius, t + timeWindowRadius + 1);

			const integer tDelta = t - pointSet.timeBegin();
			const integer tWidth = pointSet.timeEnd() - pointSet.timeBegin();
			
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
				randomAccessRange(pointSet.begin() + tDelta * trials, 
				pointSet.begin() + (tDelta + 1) * trials),
				kNearest - 1,
				kNearest, 
				0,
				&distanceArray,
				randomAccessRange(constantIterator(infinity<real>()), trials),
				0,
				entropyAlgorithm.normBijection());

			// After we have found the distances, we simply evaluate
			// the generic entropy estimator over the samples of
			// the current time instant.

			integer acceptedSamples = 0;
			real estimate = 0;
			for (integer i = 0;i < trials;++i)
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
				estimateSet[t] = entropyAlgorithm.finishEstimate(
					estimate, dimension, kNearest, 
					acceptedSamples, tWidth * trials);
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
		typename EntropyAlgorithm, 
		typename Real_OutputIterator>
	integer temporalGenericEntropy(
		const SignalPtr& signal,
		const EntropyAlgorithm& entropyAlgorithm,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		return Tim::temporalGenericEntropy(
			forwardRange(constantIterator(signal)),
			entropyAlgorithm,
			timeWindowRadius,
			result,
			kNearest);
	}

	// Generic entropy
	// ---------------

	template <
		typename SignalPtr_Iterator,
		typename EntropyAlgorithm>
	real genericEntropy(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const EntropyAlgorithm& entropyAlgorithm,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);

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

		Array<real, 2> distanceArray(1, estimateSamples);

		// Find the distance to the k:th nearest neighbor for all points.

		searchAllNeighbors(
			pointSet.kdTree(),
			randomAccessRange(pointSet.begin(), pointSet.end()),
			kNearest - 1,
			kNearest, 
			0,
			&distanceArray,
			randomAccessRange(constantIterator(infinity<real>()), estimateSamples),
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
				estimate, dimension, kNearest, 
				acceptedSamples, estimateSamples);
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

	template <typename EntropyAlgorithm>
	real genericEntropy(
		const SignalPtr& signal,
		const EntropyAlgorithm& entropyAlgorithm,
		integer kNearest)
	{
		return Tim::genericEntropy(
			forwardRange(constantIterator(signal)),
			entropyAlgorithm,
			kNearest);
	}

}

#endif
