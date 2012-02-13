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

#include <pastel/math/normbijections.h>

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/iterator_range.h>
#include <pastel/sys/eps.h>
#include <pastel/sys/copy_n.h>

#include <numeric>
#include <iterator>

#include <vector>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator,
		typename LocalEstimator>
	real entropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet,
		const ForwardIterator_Range<Integer_Iterator>& lagSet,
		integer kNearest,
		const LocalEstimator& localEstimator)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(lagSet.size(), ==, signalSet.height());

		typedef typename LocalEstimator::Instance Estimator;
		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

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

		const integer n = samples * trials;
		const integer marginals = rangeSet.size();

		// Construct point sets

		SignalPointSet jointPointSet(
			range(jointSignalSet.begin(), jointSignalSet.end()));
	
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
				Pastel::range(jointSignalSet.begin(), jointSignalSet.end()), 
				offsetSet[range[0]], offsetSet[range[1]])));
			weightSet.push_back(range[2]);
			++iter;
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_NormBijection<real> normBijection;

		// Start estimation.

		Array<real> distanceArray(1, n);

		searchAllNeighbors(
			jointPointSet.kdTree(),
			range(jointPointSet.begin(), jointPointSet.end()),
			kNearest - 1,
			kNearest, 
			(Array<Point_ConstIterator>*)0,
			&distanceArray,
			constantRange(infinity<real>(), n),
			0,
			normBijection);

		const real signalWeightSum = 
			std::accumulate(weightSet.begin(), weightSet.end(), (real)0);

		Estimator estimator(kNearest, n);

		std::vector<integer> countSet(n, 0);

		real estimate = estimator.localJointEstimate();
		for (integer i = 0;i < marginals;++i)
		{
			// Note: the maximum norm bijection values coincide 
			// with the norm values, so no need to convert.
			countAllNeighbors(
				pointSet[i]->kdTree(),
				range(pointSet[i]->begin(), pointSet[i]->end()),
				range(distanceArray.begin(), n),
				countSet.begin(),
				8,
				normBijection);
			
			integer acceptedSamples = 0;
			real signalEstimate = 0;
#pragma omp parallel for reduction(+ : signalEstimate, acceptedSamples)
			for (integer j = 0;j < n;++j)
			{
				const integer k = countSet[j];
				// A neighbor count of zero can happen when the distance
				// to the k:th neighbor is zero because of using an
				// open search ball. These points are ignored.
				if (k > 0)
				{
					signalEstimate += estimator.localMarginalEstimate(k);
					++acceptedSamples;
				}
			}
			if (acceptedSamples > 0)
			{
				signalEstimate /= acceptedSamples;
			}

			estimate -= signalEstimate * weightSet[i];
		}

		return estimate;
	}

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator>
	real entropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet,
		const ForwardIterator_Range<Integer_Iterator>& lagSet,
		integer kNearest)
	{
		return Tim::entropyCombination(
			signalSet, rangeSet,
			lagSet, kNearest,
			Log_LocalEstimator());
	}

	template <
		typename Integer3_Iterator,
		typename Real_OutputIterator>
	real entropyCombination(
		const Array<SignalPtr>& signalSet,
		const ForwardIterator_Range<Integer3_Iterator>& rangeSet)
	{
		return Tim::entropyCombination(
			signalSet,
			rangeSet,
			constantRange(0, signalSet.height()));
	}

}


#endif
