#ifndef TIM_ENTROPY_COMBINATION_T_HPP
#define TIM_ENTROPY_COMBINATION_T_HPP

#include "tim/core/entropy_combination_t.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <pastel/math/normbijection/maximum_normbijection.h>

#include <pastel/sys/iterator/constant_iterator.h>
#include <pastel/sys/iterator/counting_iterator.h>
#include <pastel/sys/range.h>
#include <pastel/sys/math/eps.h>
#include <pastel/sys/sequence/copy_n.h>
#include <pastel/sys/indicator/predicate_indicator.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

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
	Signal temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const boost::iterator_range<Integer_Iterator>& lagSet,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
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
			return Signal(0, 1);
		}

		// Check that the trials of signals have the 
		// same dimension.

		integer signals = signalSet.height();
		for (integer i = 0;i < signals;++i)
		{
			PENSURE(equalDimension(
				range(countingIterator(signalSet.cRowBegin(i)), countingIterator(signalSet.cRowEnd(i)))));
		}

		integer trials = signalSet.width();

		// Find out the shared time interval that
		// the signals share.

		Integer2 sharedTime = sharedTimeInterval(
			countingRange(signalSet.cbegin(), signalSet.cend()),
			lagSet);

		integer estimateBegin = sharedTime[0];
		integer estimateEnd = sharedTime[1];
		integer estimates = estimateEnd - estimateBegin;

		if (estimates == 0)
		{
			// The signals do not share any time interval:
			// return an empty signal.
			return Signal(0, 1);
		}
		
		// Construct the joint signal.

		std::vector<Signal> jointSignalSet;
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

		Maximum_NormBijection<real> normBijection;
		integer missingValues = 0;

		// This is where the estimates are stored at.
		
		Signal result(estimates, 1, estimateBegin);
		
		std::vector<Integer3> copyRangeSet(
			rangeSet.begin(), rangeSet.end());

		// Copy the filter and replicate
		// the values to each trial.

		integer filterWidth = filter.size();
		integer filterRadius = filterWidth / 2;
		integer maxLocalFilterWidth = 
			std::min(filterWidth, estimates);

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

#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSet jointPointSet(
			countingRange(jointSignalSet.begin(), jointSignalSet.end()));

		std::vector<SignalPointSet> pointSet;
		pointSet.reserve(marginals);

		real signalWeightSum = 0;

		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];
			
			pointSet.emplace_back(
				Pastel::countingRange(jointSignalSet.begin(), jointSignalSet.end()),
				offsetSet[range[0]], offsetSet[range[1]]);

			signalWeightSum += range[2];
		}


		Array<real> distanceArray(Vector2i(1, maxLocalFilterWidth * trials), infinity<real>());
		
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

					Vector<real> queryPoint(
						ofDimension(jointPointSet.dimension()),
						withAliasing((real*)(query->point())));

					distanceArray(i - searchBegin) =
						searchNearest(
							kdTreeNearestSet(jointPointSet.kdTree()),
							queryPoint,
							PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
							PASTEL_TAG(normBijection), normBijection,
							PASTEL_TAG(kNearest), kNearest
						).first;
				}
			};

			tbb::parallel_for(
				Block(searchBegin, searchEnd),
				search);

			real estimate = 0;
			for (integer i = 0;i < marginals;++i)
			{
				pointSet[i].setTimeWindow(
					t - timeWindowRadius, 
					t + timeWindowRadius + 1);

				real signalEstimate = 0;
				real weightSum = 0;
				integer filterOffset = tFilterOffset * trials;

				for (integer j = 0;j < windowSamples;++j)
				{
					Point_ConstIterator query =
						*(pointSet[i].begin() + searchBegin + j);

					Vector<real> queryPoint(
						ofDimension(pointSet[i].dimension()),
						withAliasing((real*)(query->point())));

					integer k = 0;
					searchNearest(
						kdTreeNearestSet(pointSet[i].kdTree()),
						queryPoint,
						PASTEL_TAG(kNearest), (integer)Infinity(),
						PASTEL_TAG(report), [&](auto, auto) {++k;},
						PASTEL_TAG(normBijection), normBijection,
						PASTEL_TAG(maxDistance2), distanceArray(j)
					);

					// Note: k = 0 is possible: a range count of zero 
					// can happen when the distance to the k:th neighbor is 
					// zero because of using an open search ball. 
									
					// These singular cases must be taken into account and
					// gracefully ignored, as is done here.

					if (k > 0)
					{
						real weight = copyFilter[j + filterOffset];
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
					estimate = (real)Nan();
					++missingValues;
					
					// Skip to the next time instant.
					break;
				}
			}

			const integer estimateSamples = tWidth * trials;

			estimate += digamma<real>(kNearest);
			estimate += (signalWeightSum - 1) * digamma<real>(estimateSamples);

			result.data()(t - estimateBegin) = estimate;
		}
		}

		// Reconstruct the NaN's in the estimates.

		reconstruct(
			range(result.data().begin(), estimates));

		return result;
	}

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	Signal temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		const boost::iterator_range<Integer_Iterator>& lagSet,
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
	Signal temporalEntropyCombination(
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
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
