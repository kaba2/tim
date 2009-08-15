#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_HPP
#define TIM_MUTUAL_INFORMATION_KRASKOW_HPP

#include "tim/core/mutual_information_kraskow.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/slidingmidpoint_splitrule_pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_neighbors_pointkdtree.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/countingiterator.h>
#include <pastel/sys/forwardrange.h>

#include <deque>

namespace Tim
{

	template <
		typename Signal_A_Iterator,
		typename Signal_B_Iterator,
		typename Real_OutputIterator>
	void mutualInformation(
		const ForwardRange<Signal_A_Iterator>& aSignalSet,
		const ForwardRange<Signal_B_Iterator>& bSignalSet,
		Real_OutputIterator result,
		integer bLag,
		integer sigma,
		integer kNearest)
	{
		ENSURE_OP(bLag, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		PENSURE_OP(aSignalSet.size(), ==, bSignalSet.size());
		PENSURE(equalDimension(aSignalSet));
		PENSURE(equalDimension(bSignalSet));

		if (aSignalSet.empty())
		{
			return;
		}

		aSignalSet.updateCache();
		bSignalSet.updateCache();

		const integer trials = aSignalSet.size();

		enum
		{
			Signals = 2
		};

		// Determine start index for each signal.

		Vector<integer, Signals> signalLag(0, bLag);
		const integer signalMaxLag = max(signalLag);
		Vector<integer, Signals> signalStartIndex(signalMaxLag - signalLag);

		// Determine the number of samples.

		Vector<integer, Signals> signalSamples(
			minSamples(aSignalSet) - signalLag[0],
			minSamples(bSignalSet) - signalLag[1]);

		const integer samples = min(signalSamples);

		if (samples <= 0)
		{
			return;
		}

		if (sigma < 0)
		{
			sigma = samples;
		}

		Infinity_NormBijection<real> normBijection;

		// Form the joint signal.

		std::vector<SignalPtr> jointSignalSet;
		jointSignalSet.reserve(trials);

		merge(aSignalSet, bSignalSet, 
			std::back_inserter(jointSignalSet), bLag);

		// Compute SignalPointSets.

		std::vector<real> estimateSet(samples);
#pragma omp parallel
		{
		SignalPointSetPtr jointPointSet(
			new SignalPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end())));

		Tuple<SignalPointSetPtr, 2> pointSet(
			SignalPointSetPtr(new SignalPointSet(aSignalSet)),
			SignalPointSetPtr(new SignalPointSet(bSignalSet)));

		Array<real, 2> distanceArray(1, trials);
		std::vector<integer> countSet(trials, 0);

#pragma omp for
		for (integer t = 0;t < samples;++t)
		{
			const integer tLeft = std::max(t - sigma, 0);
			const integer tRight = std::min(t + sigma + 1, samples);
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

			for (integer i = 0;i < Signals;++i)
			{
				pointSet[i]->setTimeWindow(tLeft, tRight);

				countAllNeighbors(
					pointSet[i]->kdTree(),
					randomAccessRange(pointSet[i]->begin() + tDelta * trials, 
					pointSet[i]->begin() + (tDelta + 1) * trials),
					randomAccessRange(distanceArray.begin(), trials),
					normBijection,
					countSet.begin());

//#pragma omp parallel for reduction(+ : estimate)
				for (integer j = 0;j < trials;++j)
				{
					estimate -= digamma<real>(countSet[j]);
				}
			}

			estimate /= trials;
			estimate += digamma<real>(kNearest);
			estimate += digamma<real>(tWidth * trials);

			estimateSet[t] = estimate;
		}
		}

		std::copy(estimateSet.begin(), estimateSet.end(), result);
	}

	template <
		typename Real_OutputIterator>
	void mutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		Real_OutputIterator result,
		integer bLag,
		integer sigma,
		integer kNearest)
	{
		SignalPtr signalSet[2] = {aSignal, bSignal};

		Tim::mutualInformation(
			forwardRange(constantIterator(aSignal)),
			forwardRange(constantIterator(bSignal)),
			result,
			bLag,
			sigma,
			kNearest);
	}

}

#endif
