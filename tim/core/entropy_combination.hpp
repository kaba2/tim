#ifndef TIM_ENTROPY_COMBINATION_HPP
#define TIM_ENTROPY_COMBINATION_HPP

#include "tim/core/entropy_combination.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/eps.h>

#include <numeric>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator,
		typename Integer_Iterator>
	void temporalEntropyCombination(
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
			return;
		}

		rangeSet.updateCache();

		const integer signals = signalSet.height();
		for (integer i = 0;i < signals;++i)
		{
			PENSURE(equalDimension(
				forwardRange(signalSet.rowBegin(i), signalSet.rowEnd(i))));
		}

		const integer trials = signalSet.width();

		// Construct the joint signal.

		std::vector<SignalPtr> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, std::back_inserter(jointSignalSet), lagSet);

		const integer samples = jointSignalSet.front()->samples();
		if (samples == 0)
		{
			return;
		}

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

		std::vector<real> estimateSet(samples);
#pragma omp parallel
		{
		// Compute SignalPointSets.

		SignalPointSetPtr jointPointSet(new SignalPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end())));

		std::vector<Integer3> copyRangeSet(
			rangeSet.begin(), rangeSet.end());

		std::vector<integer> marginalOffsetSet;
		marginalOffsetSet.reserve(marginals);

		std::vector<SignalPointSetPtr> pointSet;
		pointSet.reserve(marginals);

		real sumWeight = 0;
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = copyRangeSet[i];

			pointSet.push_back(
				SignalPointSetPtr(new SignalPointSet(
				forwardRange(jointSignalSet.begin(), jointSignalSet.end()), false,
				offsetSet[range[0]], offsetSet[range[1]])));

			sumWeight += range[2];
		}

		Array<real, 2> distanceArray(1, trials);
		std::vector<integer> countSet(trials, 0);

#pragma omp for
		for (integer t = 0;t < samples;++t)
		{
			jointPointSet->setTimeWindow(
				t - timeWindowRadius, 
				t + timeWindowRadius + 1);
			
			const integer tDelta = t - jointPointSet->timeBegin();
			const integer tWidth = jointPointSet->timeEnd() - jointPointSet->timeBegin();

			searchAllNeighbors(
				jointPointSet->kdTree(),
				DepthFirst_SearchAlgorithm_PointKdTree(),
				randomAccessRange(
				jointPointSet->begin() + tDelta * trials, 
				jointPointSet->begin() + (tDelta + 1) * trials),
				kNearest - 1,
				kNearest, 
				randomAccessRange(constantIterator(infinity<real>()), trials),
				0,
				normBijection,
				0,
				&distanceArray);

			for (integer j = 0;j < trials;++j)
			{
				distanceArray(j) = nextSmaller(distanceArray(j));
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
					pointSet[i]->begin() + tDelta * trials, 
					pointSet[i]->begin() + (tDelta + 1) * trials),
					randomAccessRange(distanceArray.begin(), trials),
					normBijection,
					countSet.begin());
				
				integer accepted = 0;
				real signalEstimate = 0;
				for (integer j = 0;j < trials;++j)
				{
					const integer k = countSet[j];
					if (k > 0)
					{
						signalEstimate += digamma<real>(k);
						++accepted;
					}
				}
				signalEstimate /= trials;

				estimate -= signalEstimate * copyRangeSet[i][2];
			}

			const integer estimateSamples = tWidth * trials;

			estimate += digamma<real>(kNearest);
			estimate += (sumWeight - 1) * digamma<real>(estimateSamples);

			estimateSet[t] = estimate;
		}
		}

		std::copy(estimateSet.begin(), estimateSet.end(), result);
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	void temporalEntropyCombination(
		const SignalPtr& signal,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		Tim::temporalEntropyCombination(
			forwardRange(constantIterator(signal)),
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
			DepthFirst_SearchAlgorithm_PointKdTree(),
			randomAccessRange(jointPointSet.begin(), jointPointSet.end()),
			kNearest - 1,
			kNearest, 
			randomAccessRange(constantIterator(infinity<real>()), estimateSamples),
			0,
			normBijection,
			0,
			&distanceArray);

#pragma omp parallel for
		for (integer j = 0;j < estimateSamples;++j)
		{
			distanceArray(j) = nextSmaller(distanceArray(j));
		}

		ENSURE_OP(*std::max_element(distanceArray.begin(), distanceArray.end()), !=, infinity<real>());

		const real sumWeight = 
			std::accumulate(weightSet.begin(), weightSet.end(), (real)0);

		std::vector<integer> countSet(estimateSamples, 0);

		real estimate = 0;
		for (integer i = 0;i < marginals;++i)
		{
			countAllNeighbors(
				pointSet[i]->kdTree(),
				randomAccessRange(pointSet[i]->begin(), pointSet[i]->end()),
				randomAccessRange(distanceArray.begin(), estimateSamples),
				normBijection,
				countSet.begin());

			integer accepted = 0;
			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate, accepted)
			for (integer j = 0;j < estimateSamples;++j)
			{
				const integer k = countSet[j];
				if (k > 0)
				{
					signalEstimate += digamma<real>(k);
					++accepted;
				}
			}
			signalEstimate /= accepted;

			estimate -= signalEstimate * weightSet[i];
		}
		
		estimate += digamma<real>(kNearest);
		estimate += (sumWeight - 1) * digamma<real>(estimateSamples);

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
			forwardRange(constantIterator(0), signalSet.height()),
			kNearest);
	}

}

#endif
