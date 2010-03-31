#ifndef TIM_ENTROPY_COMBINATION_HPP
#define TIM_ENTROPY_COMBINATION_HPP

#include "tim/core/entropy_combination.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/search_all_neighbors_pointkdtree.h>
#include <pastel/geometry/count_all_range_pointkdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/math/normbijection.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/forwardrange.h>
#include <pastel/sys/eps.h>
#include <pastel/sys/stdext_copy_n.h>

#include <numeric>
#include <iterator>

#include <vector>

namespace Tim
{

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

		const real signalWeightSum = 
			std::accumulate(weightSet.begin(), weightSet.end(), (real)0);

		std::vector<integer> countSet(estimateSamples, 0);

		real estimate = 0;
		for (integer i = 0;i < marginals;++i)
		{
			// Note: the maximum norm bijection values coincide 
			// with the norm values, so no need to convert.
			countAllRange(
				pointSet[i]->kdTree(),
				randomAccessRange(pointSet[i]->begin(), pointSet[i]->end()),
				randomAccessRange(distanceArray.begin(), estimateSamples),
				countSet.begin());
			
			integer acceptedSamples = 0;
			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate, acceptedSamples)
			for (integer j = 0;j < estimateSamples;++j)
			{
				const integer k = countSet[j] - 1;
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
