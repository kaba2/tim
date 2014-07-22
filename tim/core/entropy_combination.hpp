#ifndef TIM_ENTROPY_COMBINATION_HPP
#define TIM_ENTROPY_COMBINATION_HPP

#include "tim/core/entropy_combination.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree.h>

#include <pastel/math/maximum_normbijection.h>

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/range.h>
#include <pastel/sys/eps.h>
#include <pastel/sys/copy_n.h>
#include <pastel/sys/predicate_indicator.h>

#include <numeric>
#include <iterator>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace Tim
{

	template <
		typename Integer3_Iterator,
		typename Integer_Iterator,
		typename LocalEstimator>
	real entropyCombination(
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		const boost::iterator_range<Integer_Iterator>& lagSet,
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

		std::vector<Signal> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, std::back_inserter(jointSignalSet), lagSet);

		const integer samples = jointSignalSet.front().samples();
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
			offsetSet.push_back(offsetSet[i - 1] + signalSet(0, i - 1).dimension());
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

		// Find the distances to the k:th nearest neighbors.

		Array<real> distanceArray(Vector2i(1, n));

		using Block = tbb::blocked_range<integer>;

		auto search = [&](const Block& block)
		{
			for (integer i = block.begin(); i < block.end(); ++i)
			{
				distanceArray(i) =
					searchNearest(
					jointPointSet.kdTree(),
					*(jointPointSet.begin() + i),
					nullOutput(),
					predicateIndicator(*(jointPointSet.begin() + i), NotEqualTo()),
					normBijection)
					.kNearest(kNearest);
			}
		};

		tbb::parallel_for(Block(0, n), search);

		//const real signalWeightSum = 
		//	std::accumulate(weightSet.begin(), weightSet.end(), (real)0);

		Estimator estimator(kNearest, n);

		real estimate = estimator.localJointEstimate();
		for (integer i = 0;i < marginals;++i)
		{
			using Block = tbb::blocked_range<integer>;
			using Pair = std::pair<real, integer>;
			
			auto compute = [&](
				const Block& block,
				const Pair& start)
			{
				real signalEstimate = start.first;
				integer acceptedSamples = start.second;
				for (integer j = block.begin();j < block.end();++j)
				{
					integer k = countNearest(
						pointSet[i]->kdTree(),
						*(pointSet[i]->begin() + j),
						allIndicator(),
						normBijection)
						.maxDistance(distanceArray(j));

					// A neighbor count of zero can happen when the distance
					// to the k:th neighbor is zero because of using an
					// open search ball. These points are ignored.
					if (k > 0)
					{
						signalEstimate += estimator.localMarginalEstimate(k);
						++acceptedSamples;
					}
				}

				return Pair(signalEstimate, acceptedSamples);
			};
			
			auto reduce = [](const Pair& left, const Pair& right)
			{
				return Pair(
					left.first + right.first, 
					left.second + right.second);
			};

			real signalEstimate = 0;
			integer acceptedSamples = 0;

			std::tie(signalEstimate, acceptedSamples) = 
				tbb::parallel_reduce(
					Block(0, n),
					Pair(0, 0),
					compute,
					reduce);

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
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet,
		const boost::iterator_range<Integer_Iterator>& lagSet,
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
		const Array<Signal>& signalSet,
		const boost::iterator_range<Integer3_Iterator>& rangeSet)
	{
		return Tim::entropyCombination(
			signalSet,
			rangeSet,
			constantRange(0, signalSet.height()));
	}

}


#endif
