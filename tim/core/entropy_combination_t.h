// Description: Temporal estimation of entropy combinations

#ifndef TIM_ENTROPY_COMBINATION_T_H
#define TIM_ENTROPY_COMBINATION_T_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/sys/range.h>
#include <pastel/sys/array/array.h>
#include <pastel/sys/math/eps.h>
#include <pastel/sys/sequence/copy_n.h>
#include <pastel/sys/indicator/predicate_indicator.h>

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <pastel/math/normbijection/maximum_normbijection.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <numeric>
#include <iterator>

#include <vector>
#include <deque>

namespace Tim
{

	//! Computes a temporal entropy combination of signals.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	odd(ranges::size(filter))

	signalSet:
	An ensemble of joint signals representing trials
	of the same experiment. Note: all the marginal signals
	share the memory with these joint signals.

	rangeSet:
	A sequence of m triples T_i = (a_i, b_i, s_i), 
	where [a_i, b_i] is an interval such that picking those 
	dimensions from the joint signal X gives the marginal 
	signal X_i. The s_i is the factor by which the differential 
	entropy of such a marginal signal is multiplied before summing
	to the end-result.

	timeWindowRadius:
	The radius of the time-window in samples to use.
	Smaller values give more temporal adaptivity,
	but increase errors.

	result:
	Temporal estimates of the entropy combination of the
	signals.

	lagSet:
	Lags to apply to each signal.

	filter:
	An array of coefficients by which to weight the results
	in the time-window. The center of the array corresponds 
	to the current time instant. The width of the array can 
	be arbitrary but must be odd. The coefficients must sum
	to a non-zero value.

	kNearest:
	The k:th nearest neighbor that is used to
	estimate entropy combination.

	Returns:
	The temporal estimates in a 1d-signal.
	*/
	template <
		ranges::forward_range Integer3_Range,
		ranges::forward_range Lag_Range,
		ranges::forward_range Filter_Range>
	SignalData temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet,
		integer timeWindowRadius,
		const Lag_Range& lagSet,
		integer kNearest,
		const Filter_Range& filter)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(ranges::size(lagSet), ==, signalSet.height());
		ENSURE(odd(ranges::size(filter)));

		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (ranges::empty(signalSet) || rangeSet.empty() || ranges::empty(filter))
		{
			// There's nothing to do.
			return SignalData();
		}

		// Check that the trials of signals have the 
		// same dimension.

		integer signals = signalSet.height();
		for (integer i = 0;i < signals;++i)
		{
			PENSURE(equalDimension(range(signalSet.cRowBegin(i), signalSet.cRowEnd(i))));
		}

		integer trials = signalSet.width();

		// Find out the shared time interval that
		// the signals share.

		Integer2 sharedTime = sharedTimeInterval(range(signalSet.cbegin(), signalSet.cend()), lagSet);

		integer estimateBegin = sharedTime[0];
		integer estimateEnd = sharedTime[1];
		integer estimates = estimateEnd - estimateBegin;

		if (estimates == 0)
		{
			// The signals do not share any time interval:
			// return an empty signal.
			return SignalData();
		}
		
		// Construct the joint signal.

		std::vector<SignalData> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, std::back_inserter(jointSignalSet), lagSet);

		integer marginals = rangeSet.size();

		// Find out the dimension ranges of the marginal
		// signals.

		std::vector<integer> offsetSet;
		offsetSet.reserve(signals + 1);
		offsetSet.push_back(0);
		for (integer i = 1;i < signals + 1;++i)
		{
			integer marginalDimension = signalSet(0, i - 1).dimension();
			offsetSet.push_back(offsetSet[i - 1] + marginalDimension);
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_Norm<dreal> norm;
		integer missingValues = 0;

		// This is where the estimates are stored at.
		
		SignalData result(estimates, 1, estimateBegin);
		
		std::vector<Integer3> copyRangeSet(
			std::begin(rangeSet), std::end(rangeSet));

		// Copy the filter and replicate
		// the values to each trial.

		integer filterWidth = ranges::size(filter);
		integer filterRadius = filterWidth / 2;
		integer maxLocalFilterWidth = 
			std::min(filterWidth, estimates);

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

#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSet jointPointSet(jointSignalSet);

		std::vector<SignalPointSet> pointSet;
		pointSet.reserve(marginals);

		dreal signalWeightSum = 0;

		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];
			
			pointSet.emplace_back(
				jointSignalSet,
				offsetSet[range[0]], offsetSet[range[1]]);

			signalWeightSum += range[2];
		}


		Array<dreal> distanceArray(Vector2i(1, maxLocalFilterWidth * trials), infinity<dreal>());
		
#pragma omp for reduction(+ : missingValues)
		for (integer t = estimateBegin;t < estimateEnd;++t)
		{
			jointPointSet.setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);
			
			integer tBegin = jointPointSet.windowBegin();
			integer tEnd = jointPointSet.windowEnd();
			integer tWidth = tEnd - tBegin;
			integer tLocalFilterBegin = std::max(t - filterRadius, tBegin) - tBegin;
			integer tLocalFilterEnd = std::min(t + filterRadius + 1, tEnd) - tBegin;
			integer tFilterDelta = tBegin - (t - filterRadius);
			integer tFilterOffset = std::max(tFilterDelta, (integer)0);

			const integer windowSamples = (tLocalFilterEnd - tLocalFilterBegin) * trials;

			using Block = tbb::blocked_range<integer>;

			integer searchBegin = tLocalFilterBegin * trials;
			integer searchEnd = tLocalFilterEnd * trials;

			auto search = [&](const Block& block)
			{
				for (integer i = block.begin(); i < block.end(); ++i)
				{
					auto query = *(jointPointSet.begin() + i);

					Vector<dreal> queryPoint(
						ofDimension(jointPointSet.dimension()),
						withAliasing((dreal*)(query->point())));

					distanceArray(i - searchBegin) =
						(dreal)searchNearest(
							kdTreeNearestSet(jointPointSet.kdTree()),
							queryPoint,
							PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
							PASTEL_TAG(norm), norm,
							PASTEL_TAG(kNearest), kNearest
						).first;
				}
			};

			tbb::parallel_for(
				Block(searchBegin, searchEnd),
				search);

			dreal estimate = 0;
			for (integer i = 0;i < marginals;++i)
			{
				pointSet[i].setTimeWindow(
					t - timeWindowRadius, 
					t + timeWindowRadius + 1);

				dreal signalEstimate = 0;
				dreal weightSum = 0;
				integer filterOffset = tFilterOffset * trials;

				for (integer j = 0;j < windowSamples;++j)
				{
					Point_ConstIterator query =
						*(pointSet[i].begin() + searchBegin + j);

					Vector<dreal> queryPoint(
						ofDimension(pointSet[i].dimension()),
						withAliasing((dreal*)(query->point())));

					integer k = 0;
					searchNearest(
						kdTreeNearestSet(pointSet[i].kdTree()),
						queryPoint,
						PASTEL_TAG(kNearest), (integer)Infinity(),
						PASTEL_TAG(report), [&](auto, auto) {++k;},
						PASTEL_TAG(norm), norm,
						PASTEL_TAG(maxDistance2), norm(distanceArray(j))
					);

					// Note: k = 0 is possible: a range count of zero 
					// can happen when the distance to the k:th neighbor is 
					// zero because of using an open search ball. 
									
					// These singular cases must be taken into account and
					// gracefully ignored, as is done here.

					if (k > 0)
					{
						dreal weight = copyFilter[j + filterOffset];
						signalEstimate += weight * digamma<dreal>(k);
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
					estimate = (dreal)Nan();
					++missingValues;
					
					// Skip to the next time instant.
					break;
				}
			}

			const integer estimateSamples = tWidth * trials;

			estimate += digamma<dreal>(kNearest);
			estimate += (signalWeightSum - 1) * digamma<dreal>(estimateSamples);

			result.data()(t - estimateBegin) = estimate;
		}
		}

		// Reconstruct the NaN's in the estimates.

		reconstruct(range(result.data().range().begin(), result.data().range().begin() + estimates));

		return result;
	}

	//! Computes a temporal entropy combination of signals.
	/*!
	This is a convenience function that calls:

	temporalEntropyCombination(
		signalSet,
		rangeSet,
		timeWindowRadius,
		lagSet,
		kNearest,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
	template <
		ranges::forward_range Integer3_Range,
		ranges::forward_range Lag_Range>
	Signal temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet,
		integer timeWindowRadius,
		const Lag_Range& lagSet,
		integer kNearest = 1)
	{
		return temporalEntropyCombination(
			signalSet,
			rangeSet,
			timeWindowRadius,
			lagSet,
			kNearest,
			constantRange((dreal)1, 1));
	}

	//! Computes a temporal entropy combination of signals.
	/*!
	This is a convenience function that calls:

	temporalEntropyCombination(
		signalSet,
		rangeSet,
		timeWindowRadius,
		constantRange(0, signalSet.height()));

	See the documentation for that function.
	*/
	template <ranges::forward_range Integer3_Range>
	Signal temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet,
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
