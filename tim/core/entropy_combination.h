// Description: Estimation of entropy combinations

#ifndef TIM_ENTROPY_COMBINATION_H
#define TIM_ENTROPY_COMBINATION_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"
#include "tim/core/localestimators.h"
#include "tim/core/signalpointset.h"
#include "tim/core/reconstruction.h"

#include <pastel/geometry/pointkdtree/pointkdtree.h>
#include <pastel/geometry/search_nearest.h>
#include <pastel/geometry/nearestset/kdtree_nearestset.h>

#include <pastel/math/normbijection/maximum_normbijection.h>

#include <pastel/sys/array/array.h>
#include <pastel/sys/range.h>
#include <pastel/sys/math/eps.h>
#include <pastel/sys/sequence/copy_n.h>
#include <pastel/sys/indicator/predicate_indicator.h>

#include <numeric>
#include <iterator>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

namespace Tim
{

	//! Computes an entropy combination of signals.
	/*!
	Preconditions:
	kNearest > 0

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

	kNearest:
	The k:th nearest neighbor that is used to
	estimate entropy combination.

	estimator:
	The local estimator to use in the algorithm.
	See localestimators.txt.

	Returns:
	An estimate of the entropy combination of the signals.
	*/
	template <
		ranges::forward_range Integer3_Range,
		ranges::forward_range Lag_Range,
		typename LocalEstimator>
	dreal entropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet,
		const Lag_Range& lagSet,
		integer kNearest,
		const LocalEstimator& localEstimator)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(ranges::size(lagSet), ==, signalSet.height());

		typedef typename LocalEstimator::Instance Estimator;
		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (ranges::empty(signalSet) || rangeSet.empty())
		{
			return 0;
		}

		// Construct the joint signal.

		integer trials = signalSet.width();

		std::vector<Signal> jointSignalSet;
		jointSignalSet.reserve(trials);
		merge(signalSet, 
			std::back_inserter(jointSignalSet), lagSet);

		integer samples = std::begin(signalSet)->samples();
		if (samples == 0)
		{
			return 0;
		}

		integer signals = signalSet.height();

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
		integer marginals = rangeSet.size();

		// Construct point sets

		SignalPointSet jointPointSet(jointSignalSet);
	
		std::vector<integer> weightSet;
		weightSet.reserve(marginals);

		auto iter = std::begin(rangeSet);
		std::vector<SignalPointSet> pointSet;
		pointSet.reserve(marginals);
		for (integer i = 0;i < marginals;++i)
		{
			const Integer3& range = *iter;
			pointSet.emplace_back(
				jointSignalSet,
				offsetSet[range[0]], offsetSet[range[1]]);
			weightSet.push_back(range[2]);
			++iter;
		}

		// It is essential that the used norm is the
		// maximum norm.

		Maximum_Norm<dreal> norm;

		// Find the distances to the k:th nearest neighbors.

		Array<dreal> distanceArray(Vector2i(1, n));

		using Block = tbb::blocked_range<integer>;

		auto search = [&](const Block& block)
		{
			for (integer i = block.begin(); i < block.end(); ++i)
			{
				auto query = *(jointPointSet.begin() + i);

				Vector<dreal> queryPoint(
					ofDimension(jointPointSet.dimension()),
					withAliasing((dreal*)(query->point())));

				distanceArray(i) =
					(dreal)searchNearest(
						kdTreeNearestSet(jointPointSet.kdTree()),
						queryPoint,
						PASTEL_TAG(accept), predicateIndicator(query, NotEqualTo()),
						PASTEL_TAG(norm), norm,
						PASTEL_TAG(kNearest), kNearest
					).first;
			}
		};

		tbb::parallel_for(Block(0, n), search);

		//const dreal signalWeightSum = 
		//	std::accumulate(weightSet.begin(), weightSet.end(), (dreal)0);

		Estimator estimator(kNearest, n);

		dreal estimate = estimator.localJointEstimate();
		for (integer i = 0;i < marginals;++i)
		{
			using Block = tbb::blocked_range<integer>;
			using Pair = std::pair<dreal, integer>;
			
			auto compute = [&](
				const Block& block,
				const Pair& start)
			{
				dreal signalEstimate = start.first;
				integer acceptedSamples = start.second;
				for (integer j = block.begin();j < block.end();++j) 
				{
					auto query = *(pointSet[i].begin() + j);

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

			dreal signalEstimate = 0;
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

	//! Computes an entropy combination of signals.
	/*!
	This is a convenience function that calls:

	entropyCombination(
		signalSet,
		rangeSet,
		constantRange(0, signalSet.height()),
		kNearest,
		Log_LocalEstimator());

	See the documentation for that function.
	*/
	template <
		ranges::forward_range Integer3_Range,
		ranges::forward_range Lag_Range>
	dreal entropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet,
		const Lag_Range& lagSet,
		integer kNearest = 1)
	{
		return Tim::entropyCombination(
			signalSet, rangeSet,
			lagSet, kNearest,
			Log_LocalEstimator());
	}

	//! Computes an entropy combination of signals.
	/*!
	This is a convenience function that calls:

	entropyCombination(
		signalSet,
		rangeSet,
		constantRange(0, signalSet.height()));

	See the documentation for that function.
	*/
	template <ranges::forward_range Integer3_Range>
	dreal entropyCombination(
		const Array<Signal>& signalSet,
		const Integer3_Range& rangeSet)
	{
		return Tim::entropyCombination(
			signalSet,
			rangeSet,
			constantRange(0, signalSet.height()));
	}

}

#endif
