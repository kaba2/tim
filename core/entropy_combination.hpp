#ifndef TIM_ENTROPY_COMBINATION_HPP
#define TIM_ENTROPY_COMBINATION_HPP

#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/forwardrange.h>

namespace Tim
{

	template <
		typename Signal_Iterator,
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	void temporalEntropyCombination(
		const ForwardRange<Signal_Iterator>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer kNearest)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		PENSURE(equalDimension(signalSet));

		if (signalSet.empty() || rangeSet.empty())
		{
			return;
		}

		signalSet.updateCache();
		rangeSet.updateCache();

		const integer trials = signalSet.size();
		const integer signals = rangeSet.size();
		const integer samples = minSamples(signalSet);

		if (samples == 0)
		{
			return;
		}

		Infinity_NormBijection<real> normBijection;

		// Compute SignalPointSets.

		std::vector<real> estimateSet(samples);
#pragma omp parallel
		{
		SignalPointSetPtr jointPointSet(
			new SignalPointSet(
			forwardRange(signalSet.begin(), signalSet.end())));

		std::vector<integer> weightSet;
		weightSet.reserve(signals);

		Integer3_Iterator iter = rangeSet.begin();
		std::vector<SignalPointSetPtr> pointSet;
		pointSet.reserve(signals);
		for (integer i = 0;i < signals;++i)
		{
			const Integer3& range = *iter;

			pointSet.push_back(
				SignalPointSetPtr(new SignalPointSet(
				signalSet, SignalPointSet_TimeWindow::StartEmpty, 
				range[0], range[1])));

			weightSet.push_back(range[2]);

			++iter;
		}

		Array<real, 2> distanceArray(1, trials);
		std::vector<integer> countSet(trials, 0);

#pragma omp for
		for (integer t = 0;t < samples;++t)
		{
			const integer tLeft = std::max(t - timeWindowRadius, 0);
			const integer tRight = std::min(t + timeWindowRadius + 1, samples);
			const integer tDelta = t - tLeft;
			const integer tWidth = tRight - tLeft;

			jointPointSet->setTimeWindow(tLeft, tRight);
			
			searchAllNeighbors(
				jointPointSet->kdTree(),
				DepthFirst_SearchAlgorithm_PointKdTree(),
				randomAccessRange(jointPointSet->begin() + tDelta * trials, 
				jointPointSet->begin() + (tDelta + 1) * trials),
				kNearest - 1,
				kNearest, 
				randomAccessRange(constantIterator(infinity<real>()), trials),
				0,
				normBijection,
				0,
				&distanceArray);

			real estimate = 0;

			for (integer i = 0;i < signals;++i)
			{
				pointSet[i]->setTimeWindow(tLeft, tRight);

				countAllNeighbors(
					pointSet[i]->kdTree(),
					randomAccessRange(pointSet[i]->begin() + tDelta * trials, 
					pointSet[i]->begin() + (tDelta + 1) * trials),
					randomAccessRange(distanceArray.begin(), trials),
					normBijection,
					countSet.begin());
				
				real signalEstimate = 0;
//#pragma omp parallel for reduction(+ : estimate)
				for (integer j = 0;j < trials;++j)
				{
					signalEstimate -= digamma<real>(countSet[j]);
				}

				estimate += signalEstimate * weightSet[i];
			}

			const integer estimateSamples = tWidth * trials;

			estimate /= estimateSamples;
			estimate += digamma<real>(kNearest);
			estimate += digamma<real>(estimateSamples);

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
		Tim::mutualInformation(
			forwardRange(constantIterator(signal)),
			rangeSet,
			timeWindowRadius,
			result,
			kNearest);
	}

	template <
		typename Signal_Iterator,
		typename Integer3_Iterator>
	real entropyCombination(
		const ForwardRange<Signal_Iterator>& signalSet,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer kNearest)
	{
		ENSURE_OP(kNearest, >, 0);

		const integer trials = signalSet.size();
		const integer samples = minSamples(signalSet);
		const integer dimension = signalSet.front()->dimension();
		const integer estimateSamples = samples * trials;
		const integer signals = rangeSet.size();

		// Construct point sets

		SignalPointSet jointPointSet(signalSet, SignalPointSet_TimeWindow::StartFull);
	
		std::vector<integer> weightSet;
		weightSet.reserve(signals);

		Integer3_Iterator iter = rangeSet.begin();
		std::vector<SignalPointSetPtr> pointSet;
		pointSet.reserve(signals);
		for (integer i = 0;i < signals;++i)
		{
			const Integer3& range = *iter;
			pointSet.push_back(
				SignalPointSetPtr(new SignalPointSet(
				signalSet, SignalPointSet_TimeWindow::StartEmpty, 
				range[0], range[1])));
			weightSet.push_back(range[2]);
			++iter;
		}

		Infinity_NormBijection<real> normBijection;
		Array<real, 2> distanceArray(1, estimateSamples);

		// Start estimation.

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

		real estimate = 0;
		std::vector<integer> countSet(estimateSamples, 0);
		for (integer i = 0;i < signals;++i)
		{
			countAllNeighbors(
				pointSet[i]->kdTree(),
				randomAccessRange(pointSet[i]->begin(), pointSet[i]->end()),
				randomAccessRange(distanceArray.begin(), estimateSamples),
				normBijection,
				countSet.begin());

			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate)
			for (integer j = 0;j < estimateSamples;++j)
			{
				signalEstimate -= digamma<real>(countSet[j]);
			}

			estimate += weightSet[i] * signalEstimate;
		}
		estimate /= estimateSamples;
		estimate += digamma<real>(kNearest);
		estimate += digamma<real>(estimateSamples);

		return estimate;
	}

	template <typename Integer3_Iterator>
	real entropyCombination(
		const SignalPtr& signal,
		const ForwardRange<Integer3_Iterator>& rangeSet,
		integer kNearest)
	{
		return Tim::entropyCombination(
			forwardRange(constantIterator(signal)),
			rangeSet,
			kNearest);
	}

}

#endif
