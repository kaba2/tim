#ifndef TIM_ENTROPY_COMBINATION_HPP
#define TIM_ENTROPY_COMBINATION_HPP

#include "tim/core/entropy_combination.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/eps.h>
#include <pastel/sys/stdext_copy_n.h>

#include <numeric>
#include <iterator>

#include <vector>
#include <deque>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator,
		typename Integer_Iterator>
	integer temporalEntropyCombination(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer kNearest)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(lagSet.size(), ==, signalSet.height());

		if (signalSet.empty() || rangeSet.empty())
		{
			// There's nothing to do.
			return 0;
		}

		rangeSet.updateCache();

		// Check that the trials of signals have the 
		// same dimension.

		const integer signals = signalSet.height();
		for (integer i = 0;i < signals;++i)
		{
			PENSURE(equalDimension(
				forwardRange(signalSet.rowBegin(i), signalSet.rowEnd(i))));
		}

		const integer trials = signalSet.width();

		// There are several goals to keep in mind when
		// reporting the results:
		// 1) We want the size of the output data to be 
		// independent of the lags.
		// 2) The time instant of the first element of the
		// output data must be 0.
		// 3) The width of the output data must be able
		// to contain all temporal estimates no matter
		// which lags are used (this requires the user
		// to choose lags so that the common time interval
		// of the signals does not lie outside the 
		// time-window of the output data).

		// The requirements 1 and 3 are satisfied by
		// choosing the width of the output data as
		// the minimum number of samples among all
		// involved signals.
	
		const integer outputWidth = minSamples(
			forwardRange(signalSet.begin(), signalSet.end()));

		// Find out the shared time interval that
		// the signals share.

		const Integer2 sharedTime = sharedTimeInterval(
			forwardRange(signalSet.begin(), signalSet.end()),
			lagSet);

		if (sharedTime[0] == sharedTime[1] ||
			sharedTime[0] >= outputWidth ||
			sharedTime[1] <= 0)
		{
			// The shared time interval is not in the
			// output range. Return all NaNs.
			std::fill_n(result, outputWidth, nan<real>());
			return 0;
		}
		
		// We are only interested in estimates in the
		// output range.

		const integer estimateBegin = std::max(sharedTime[0], 0);
		const integer estimateEnd = std::min(sharedTime[1], outputWidth);
		const integer estimates = estimateEnd - estimateBegin;
		const integer estimateOffset = estimateBegin - sharedTime[0];

		// Construct the joint signal.

		std::vector<SignalPtr> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, std::back_inserter(jointSignalSet), lagSet);

		const integer marginals = rangeSet.size();

		// Find out the dimension ranges of the marginal
		// signals.

		std::vector<integer> offsetSet;
		offsetSet.reserve(signals + 1);
		offsetSet.push_back(0);
		for (integer i = 1;i < signals + 1;++i)
		{
			offsetSet.push_back(offsetSet[i - 1] + signalSet(0, i - 1)->dimension());
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_NormBijection<real> normBijection;
		integer missingValues = 0;

		// The requirement 2 is satisfied by padding
		// the output data before and after the estimates
		// with NaNs.

		std::vector<real> estimateSet(outputWidth, nan<real>());

		const integer maxSamplesInTimeWindow = 
			std::min(2 * timeWindowRadius + 1, estimates);

		std::vector<Integer3> copyRangeSet(
			rangeSet.begin(), rangeSet.end());

		std::vector<real> copyWeightSet;
		copyWeightSet.reserve(maxSamplesInTimeWindow * trials);

		// Create a triangle filter.
		for (integer t = 0;t < maxSamplesInTimeWindow;++t)
		{
			const real x = 
				(real)std::abs(t - timeWindowRadius) / 
				(timeWindowRadius + (real)0.5);
			std::fill_n(
				std::back_inserter(copyWeightSet), trials, 1 - x);
		}

//#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSetPtr jointPointSet(new SignalPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end())));

		std::vector<SignalPointSetPtr> pointSet(marginals);

		real signalWeightSum = 0;
#pragma omp parallel for reduction(+ : signalWeightSum)
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];

			const SignalPointSetPtr marginalPointSet(
				new SignalPointSet(
				forwardRange(jointSignalSet.begin(), jointSignalSet.end()), false,
				offsetSet[range[0]], offsetSet[range[1]]));

			pointSet[i] = marginalPointSet;
				
			signalWeightSum += range[2];
		}

		Array<real, 2> distanceArray(1, maxSamplesInTimeWindow * trials, infinity<real>());
		std::deque<real> hintDistanceSet(
			maxSamplesInTimeWindow * trials, infinity<real>());
		
		std::vector<integer> countSet(maxSamplesInTimeWindow * trials, 0);

//#pragma omp for reduction(+ : missingValues)
		for (integer t = estimateOffset;t < estimateOffset + estimates;++t)
		{
			jointPointSet->setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);
			
			const integer tBegin = jointPointSet->timeBegin();
			const integer tEnd = jointPointSet->timeEnd();
			const integer tWidth = tEnd - tBegin;
			const integer tDelta = tBegin - (t - timeWindowRadius);

			const integer estimateSamples = tWidth * trials;
			const real maxRelativeError = 0;

			searchAllNeighbors(
				jointPointSet->kdTree(),
				randomAccessRange(
				jointPointSet->begin(), 
				jointPointSet->end()),
				kNearest - 1,
				kNearest, 
				0,
				&distanceArray,
				constantRange(infinity<real>(), estimateSamples),
				maxRelativeError,
				normBijection,
				randomAccessRange(hintDistanceSet.begin(), estimateSamples));
			
			if (tDelta == 0)
			{
				// Shift the hint distances to the left
				// for the next time instant.
				for (integer i = 0;i < trials;++i)
				{
					hintDistanceSet.push_back(infinity<real>());
					hintDistanceSet.pop_front();
				}
			}

			for (integer j = 0;j < estimateSamples;++j)
			{
				// This is an important part of the algorithm although it may
				// not seem like it from the first look. Because of the use of
				// the maximum norm, the k:th neighbor will be at least on the
				// surface of one of the marginal balls. Those points on the 
				// surfaces must not be counted or otherwise the result becomes
				// biased. I.e. one must use an open search ball for counting marginal
				// neighbors, and not a closed one. This is not a rare case: 
				// it happens in the marginal neighbor searching for every point. 
				// Tracing this bug took many days.

				distanceArray(j) = std::max(nextSmaller(distanceArray(j)), (real)0);
			}

			real estimate = 0;
			for (integer i = 0;i < marginals;++i)
			{
				pointSet[i]->setTimeWindow(
					t - timeWindowRadius, 
					t + timeWindowRadius + 1);

				countAllNeighbors(
					pointSet[i]->kdTree(),
					randomAccessRange(
					pointSet[i]->begin(), 
					pointSet[i]->end()),
					randomAccessRange(distanceArray.begin(), estimateSamples),
					countSet.begin(),
					normBijection);
				
				real signalEstimate = 0;
				real weightSum = 0;
				const integer weightOffset = tDelta * trials;
#pragma omp parallel for reduction(+ : signalEstimate, weightSum)
				for (integer j = 0;j < estimateSamples;++j)
				{
					const integer k = countSet[j];

					// A neighbor count of zero can happen when the distance
					// to the k:th neighbor is zero because of using an
					// open search ball. These points are ignored.
					if (k > 0)
					{
						const real weight = copyWeightSet[j + weightOffset];
						signalEstimate += weight * digamma<real>(k);
						weightSum += weight;
					}
				}
				if (weightSum != 0)
				{
					signalEstimate /= weightSum;
					estimate -= signalEstimate * copyRangeSet[i][2];
				}
				else
				{
					// The estimate is undefined, mark
					// it with NaN. This value will
					// probably be reconstructed later.
					estimate = nan<real>();
					++missingValues;
					
					// Skip to the next time instant.
					break;
				}
			}

			estimate += digamma<real>(kNearest);
			estimate += (signalWeightSum - 1) * digamma<real>(estimateSamples);

			estimateSet[t - estimateOffset + estimateBegin] = estimate;
		}
		}

		// Reconstruct the NaN's in the estimates.

		reconstruct(
			forwardRange(estimateSet.begin() + estimateBegin, estimates));

		// Copy the results to output.

		std::copy(estimateSet.begin(), estimateSet.end(), result);

		// Done. Report the number of undefined
		// temporal estimates.

		return missingValues;
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	integer temporalEntropyCombination(
		const SignalPtr& signal,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		return Tim::temporalEntropyCombination(
			constantRange(signal),
			rangeSet,
			timeWindowRadius,
			result,
			kNearest);
	}

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	real entropyCombination(
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
	
		std::vector<integer> weightSet;
		weightSet.reserve(marginals);

		Integer3_Iterator iter = rangeSet.begin();
		std::vector<SignalPointSetPtr> pointSet;
		pointSet.reserve(marginals);
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = *iter;
			pointSet.push_back(
				SignalPointSetPtr(new SignalPointSet(
				forwardRange(jointSignalSet.begin(), jointSignalSet.end()), 
				true,
				offsetSet[range[0]], offsetSet[range[1]])));
			weightSet.push_back(range[2]);
			++iter;
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_NormBijection<real> normBijection;

		// Start estimation.

		Array<real, 2> distanceArray(1, estimateSamples);

		searchAllNeighbors(
			jointPointSet.kdTree(),
			randomAccessRange(jointPointSet.begin(), jointPointSet.end()),
			kNearest - 1,
			kNearest, 
			0,
			&distanceArray,
			constantRange(infinity<real>(), estimateSamples),
			0,
			normBijection);

#pragma omp parallel for
		for (integer j = 0;j < estimateSamples;++j)
		{
			// This is an important part of the algorithm although it may
			// not seem like it from the first look. Because of the use of
			// the maximum norm, the k:th neighbor will be at least on the
			// surface of one of the marginal balls. Those points on the 
			// surfaces must not be counted or otherwise the result becomes
			// biased. I.e. one must use an open search ball for counting marginal
			// neighbors, and not a closed one. This is not a rare case: 
			// it happens in the marginal neighbor searching for every point. 
			// Tracing this bug took many days.

			distanceArray(j) = std::max(nextSmaller(distanceArray(j)), (real)0);
		}

		const real signalWeightSum = 
			std::accumulate(weightSet.begin(), weightSet.end(), (real)0);

		std::vector<integer> countSet(estimateSamples, 0);

		real estimate = 0;
		for (integer i = 0;i < marginals;++i)
		{
			countAllNeighbors(
				pointSet[i]->kdTree(),
				randomAccessRange(pointSet[i]->begin(), pointSet[i]->end()),
				randomAccessRange(distanceArray.begin(), estimateSamples),
				countSet.begin(),
				normBijection);

			integer acceptedSamples = 0;
			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate, acceptedSamples)
			for (integer j = 0;j < estimateSamples;++j)
			{
				const integer k = countSet[j];
				// A neighbor count of zero can happen when the distance
				// to the k:th neighbor is zero because of using an
				// open search ball. These points are ignored.
				if (k > 0)
				{
					signalEstimate += digamma<real>(k);
					++acceptedSamples;
				}
			}
			if (acceptedSamples > 0)
			{
				signalEstimate /= acceptedSamples;
			}

			estimate -= signalEstimate * weightSet[i];
		}
		
		estimate += digamma<real>(kNearest);
		estimate += (signalWeightSum - 1) * digamma<real>(estimateSamples);

		return estimate;
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	real entropyCombination(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet)
	{
		return Tim::entropyCombination(
			signalSet,
			rangeSet,
			constantRange(0, signalSet.height()),
			kNearest);
	}

}

#endif
