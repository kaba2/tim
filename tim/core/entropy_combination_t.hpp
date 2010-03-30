#ifndef TIM_ENTROPY_COMBINATION_T_HPP
#define TIM_ENTROPY_COMBINATION_T_HPP

#include "tim/core/entropy_combination_t.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_range_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>
#include <pastel/geometry/search_depth_first_pointkdtree.h>

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
		typename Integer_Iterator,
		typename Real_Filter_Iterator>
	integer temporalEntropyCombination(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		const ForwardRange<Integer_Iterator>& lagSet,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(lagSet.size(), ==, signalSet.height());
		ENSURE(odd(filter.size()));

		if (signalSet.empty() || rangeSet.empty() || filter.empty())
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

		log() << "Trials = " << trials << logNewLine;

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

		std::vector<Integer3> copyRangeSet(
			rangeSet.begin(), rangeSet.end());

		// Copy the filter and replicate
		// the values to each trial.

		const integer filterWidth = filter.size();
		const integer filterRadius = filterWidth / 2;
		const integer maxLocalFilterWidth = 
			std::min(filterWidth, estimates);

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

		/*
		// Create a triangle filter.
		for (integer t = 0;t < filterWidth;++t)
		{
			const real x = 
				(real)std::abs(t - timeWindowRadius) / 
				(timeWindowRadius + (real)0.5);
			std::fill_n(
				std::back_inserter(copyFilter), trials, 1 - x);
		}
		*/

// Parallelization version 1
#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSetPtr jointPointSet(new SignalPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end())));

		std::vector<SignalPointSetPtr> pointSet(marginals);

		real signalWeightSum = 0;
// Parallelization version 2
//#pragma omp parallel for reduction(+ : signalWeightSum)
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

		Array<real, 2> distanceArray(1, maxLocalFilterWidth * trials, infinity<real>());
		
		std::vector<integer> countSet(maxLocalFilterWidth * trials, 0);

// Parallelization version 1
#pragma omp for reduction(+ : missingValues)
		for (integer t = estimateOffset;t < estimateOffset + estimates;++t)
		{
			jointPointSet->setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);
			
			const integer tBegin = jointPointSet->timeBegin();
			const integer tEnd = jointPointSet->timeEnd();
			const integer tWidth = tEnd - tBegin;
			const integer tLocalFilterBegin = std::max(t - filterRadius, tBegin) - tBegin;
			const integer tLocalFilterEnd = std::min(t + filterRadius + 1, tEnd) - tBegin;
			const integer tFilterDelta = tBegin - (t - filterRadius);
			const integer tFilterOffset = std::max(tFilterDelta, 0);
			const integer windowSamples = (tLocalFilterEnd - tLocalFilterBegin) * trials;
			const real maxRelativeError = 0;

			searchAllNeighbors(
				jointPointSet->kdTree(),
				randomAccessRange(
				jointPointSet->begin() + tLocalFilterBegin * trials, 
				jointPointSet->begin() + tLocalFilterEnd * trials),
				kNearest - 1,
				kNearest, 
				0,
				&distanceArray,
				constantRange(infinity<real>(), windowSamples),
				maxRelativeError,
				normBijection,
				DepthFirst_SearchAlgorithm_PointKdTree());
			
			real estimate = 0;
			for (integer i = 0;i < marginals;++i)
			{
				pointSet[i]->setTimeWindow(
					t - timeWindowRadius, 
					t + timeWindowRadius + 1);

				// Note: the maximum norm bijection values coincide 
				// with the norm values, so no need to convert.
				countAllRange(
					pointSet[i]->kdTree(),
					randomAccessRange(
					pointSet[i]->begin() + tLocalFilterBegin * trials, 
					pointSet[i]->begin() + tLocalFilterEnd * trials),
					randomAccessRange(distanceArray.begin(), windowSamples),
					countSet.begin());
				
				real signalEstimate = 0;
				real weightSum = 0;
				const integer filterOffset = tFilterOffset * trials;
				/*
				printf("Filter:\n");
				for (integer j = 0;j < windowSamples;++j)
				{
					printf("%f ", copyFilter[j + filterOffset]);
				}
				printf("\n");
				*/

// Parallelization version 2
//#pragma omp parallel for reduction(+ : signalEstimate, weightSum)
				for (integer j = 0;j < windowSamples;++j)
				{
					const integer k = countSet[j] - 1;

					// A neighbor count of zero can happen when the distance
					// to the k:th neighbor is zero because of using an
					// open search ball. These points are ignored.
					if (k > 0)
					{
						const real weight = copyFilter[j + filterOffset];
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

			const integer estimateSamples = tWidth * trials;

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
		return temporalEntropyCombination(
			signalSet,
			rangeSet,
			timeWindowRadius,
			result,
			lagSet,
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	integer temporalEntropyCombination(
		const Array<SignalPtr, 2>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result)
	{
		return temporalEntropyCombination(
			signalSet,
			rangeSet,
			timeWindowRadius,
			result,
			constantRange(0, signalSet.height()));
	}
}

#endif
