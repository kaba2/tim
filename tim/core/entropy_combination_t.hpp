#ifndef TIM_ENTROPY_COMBINATION_T_HPP
#define TIM_ENTROPY_COMBINATION_T_HPP

#include "tim/core/entropy_combination_t.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/math/normbijections.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/iterator_range.h>
#include <pastel/sys/eps.h>
#include <pastel/sys/copy_n.h>

#include <numeric>
#include <iterator>

#include <vector>
#include <deque>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator,
		typename Real_Filter_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const ForwardIterator_Range<Integer_Iterator>& lagSet,
		integer kNearest,
		const ForwardIterator_Range<Real_Filter_Iterator>& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(lagSet.size(), ==, signalSet.height());
		ENSURE(odd(filter.size()));

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (signalSet.empty() || rangeSet.empty() || filter.empty())
		{
			// There's nothing to do.
			return SignalPtr(new Signal(0, 1));
		}

		rangeSet.updateCache();

		// Check that the trials of signals have the 
		// same dimension.

		const integer signals = signalSet.height();
		for (integer i = 0;i < signals;++i)
		{
			PENSURE(equalDimension(
				range(signalSet.rowBegin(i), signalSet.rowEnd(i))));
		}

		const integer trials = signalSet.width();

		// Find out the shared time interval that
		// the signals share.

		const Integer2 sharedTime = sharedTimeInterval(
			range(signalSet.begin(), signalSet.end()),
			lagSet);

		const integer estimateBegin = sharedTime[0];
		const integer estimateEnd = sharedTime[1];
		const integer estimates = estimateEnd - estimateBegin;

		if (estimates == 0)
		{
			// The signals do not share any time interval:
			// return an empty signal.
			return SignalPtr(new Signal(0, 1));
		}
		
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
			const integer marginalDimension = signalSet(0, i - 1)->dimension();
			offsetSet.push_back(offsetSet[i - 1] + marginalDimension);
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_NormBijection<real> normBijection;
		integer missingValues = 0;

		// This is where the estimates are stored at.
		
		SignalPtr result(new Signal(estimates, 1, estimateBegin));
		
		const std::vector<Integer3> copyRangeSet(
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

// Parallelization version 1
#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSetPtr jointPointSet(new SignalPointSet(
			range(jointSignalSet.begin(), jointSignalSet.end())));

		std::vector<SignalPointSetPtr> pointSet(marginals);

		real signalWeightSum = 0;
// Parallelization version 2
//#pragma omp parallel for reduction(+ : signalWeightSum)
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];

			const SignalPointSetPtr marginalPointSet(
				new SignalPointSet(
				Pastel::range(jointSignalSet.begin(), jointSignalSet.end()),
				offsetSet[range[0]], offsetSet[range[1]]));

			pointSet[i] = marginalPointSet;
				
			signalWeightSum += range[2];
		}

		Array<real> distanceArray(1, maxLocalFilterWidth * trials, infinity<real>());
		
		std::vector<integer> countSet(maxLocalFilterWidth * trials, 0);

// Parallelization version 1
#pragma omp for reduction(+ : missingValues)
		for (integer t = estimateBegin;t < estimateEnd;++t)
		{
			jointPointSet->setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);
			
			const integer tBegin = jointPointSet->windowBegin();
			const integer tEnd = jointPointSet->windowEnd();
			const integer tWidth = tEnd - tBegin;
			const integer tLocalFilterBegin = std::max(t - filterRadius, tBegin) - tBegin;
			const integer tLocalFilterEnd = std::min(t + filterRadius + 1, tEnd) - tBegin;
			const integer tFilterDelta = tBegin - (t - filterRadius);
			const integer tFilterOffset = std::max(tFilterDelta, 0);
			const integer windowSamples = (tLocalFilterEnd - tLocalFilterBegin) * trials;
			const real maxRelativeError = 0;

			searchAllNeighbors(
				jointPointSet->kdTree(),
				range(
				jointPointSet->begin() + tLocalFilterBegin * trials, 
				jointPointSet->begin() + tLocalFilterEnd * trials),
				kNearest - 1,
				kNearest, 
				(Array<Point_ConstIterator>*)0,
				&distanceArray,
				constantRange(infinity<real>(), windowSamples),
				maxRelativeError,
				normBijection);
			
			real estimate = 0;
			for (integer i = 0;i < marginals;++i)
			{
				pointSet[i]->setTimeWindow(
					t - timeWindowRadius, 
					t + timeWindowRadius + 1);

				// Note: the maximum norm bijection values coincide 
				// with the norm values, so no need to convert.
				countAllNeighbors(
					pointSet[i]->kdTree(),
					range(
					pointSet[i]->begin() + tLocalFilterBegin * trials, 
					pointSet[i]->begin() + tLocalFilterEnd * trials),
					range(distanceArray.begin(), windowSamples),
					countSet.begin(),
					8,
					normBijection);

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
					const integer k = countSet[j];

					// Note: k = 0 is possible: a range count of zero 
					// can happen when the distance to the k:th neighbor is 
					// zero because of using an open search ball. 
									
					// These singular cases must be taken into account and
					// gracefully ignored, as is done here.

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

			result->data()(t - estimateBegin) = estimate;
		}
		}

		// Reconstruct the NaN's in the estimates.

		reconstruct(
			range(result->data().begin(), estimates));

		return result;
	}

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const ForwardIterator_Range<Integer_Iterator>& lagSet,
		integer kNearest)
	{
		return temporalEntropyCombination(
			signalSet,
			rangeSet,
			timeWindowRadius,
			lagSet,
			kNearest,
			constantRange((real)1, 1));
	}

	template <typename Integer3_Iterator>
	SignalPtr temporalEntropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius)
	{
		return temporalEntropyCombination(
			signalSet,
			rangeSet,
			timeWindowRadius,
			constantRange(0, signalSet.height()));
	}
}

#endif
