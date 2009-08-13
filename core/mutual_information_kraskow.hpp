#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_HPP
#define TIM_MUTUAL_INFORMATION_KRASKOW_HPP

#include "tim/core/mutual_information_kraskow.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/pointkdtree.h>
#include <pastel/geometry/pointkdtree_refine.h>
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

		const integer trials = aSignalSet.size();

		const integer samples = std::min(
			minSamples(aSignalSet),
			minSamples(bSignalSet));

		if (sigma < 0)
		{
			sigma = samples;
		}

		const integer sigmaSamples = sigma * 2 + 1;

		const integer aDimension = aSignalSet.front()->dimension();
		const integer bDimension = bSignalSet.front()->dimension();
		const integer jointDimension = aDimension + bDimension;

		const integer maxPointsPerNode = 16;
		SlidingMidpoint2_SplitRule_PointKdTree splitRule;

		const Infinity_NormBijection<real> normBijection;
		Array<real, 2> distanceArray(1, samples);

		typedef PointKdTree<real, Dynamic, Pointer_ObjectPolicy_PointKdTree<real, Dynamic> > KdTree;
		typedef typename KdTree::ConstObjectIterator ConstObjectIterator;
		typedef typename KdTree::Object Object;

		const integer initialSamples = std::max(sigma + 1, samples);

		// Form the joint signal.

		std::vector<SignalPtr> jointSignalSet;
		jointSignalSet.reserve(trials);

		merge(aSignalSet, bSignalSet, 
			std::back_inserter(jointSignalSet), bLag);

		// Compute kd-tree for joint samples
		// ---------------------------------

		std::vector<PointD> jointPointSet;
		constructPointSet(
			forwardRange(jointSignalSet.begin(), jointSignalSet.end()),
			0, samples,
			0, jointDimension,
			jointPointSet);

		KdTree jointKdTree(ofDimension(jointDimension), maxPointsPerNode);

		jointKdTree.insert(
			countingIterator(&jointPointSet.front()), 
			countingIterator(&jointPointSet.front() + jointPointSet.size()));

		// Compute a fine subdivision for the points.

		jointKdTree.refine(splitRule);

		// Remove all objects but leave subdivision intact.

		jointKdTree.eraseObjects();

		// Insert the first time-window.

		std::deque<ConstObjectIterator> jointIteratorSet;
		for (integer i = 0;i < initialSamples;++i)
		{
			jointIteratorSet.push_back(jointKdTree.insert(&jointPointSet[i]));
		}

		// Compute kd-tree for marginal A samples
		// --------------------------------------

		std::vector<PointD> aPointSet;
		constructPointSet(
			aSignalSet,
			0, samples,
			0, aDimension,
			aPointSet);

		KdTree aKdTree(ofDimension(aDimension), maxPointsPerNode);

		aKdTree.insert(
			countingIterator(&aPointSet.front()), 
			countingIterator(&aPointSet.front() + aPointSet.size()));

		// Compute a fine subdivision for the points.

		aKdTree.refine(splitRule);

		// Remove all objects but leave subdivision intact.

		aKdTree.eraseObjects();

		// Insert the first time-window.
		
		std::deque<ConstObjectIterator> aIteratorSet;
		for (integer i = 0;i < initialSamples;++i)
		{
			aIteratorSet.push_back(aKdTree.insert(&aPointSet[i]));
		}

		// Compute kd-tree for marginal B samples
		// --------------------------------------

		std::vector<PointD> bPointSet;
		constructPointSet(
			bSignalSet, 
			0, samples,
			0, bDimension,
			bPointSet);

		KdTree bKdTree(ofDimension(bDimension), maxPointsPerNode);

		bKdTree.insert(
			countingIterator(&bPointSet.front()), 
			countingIterator(&bPointSet.front() + bPointSet.size()));

		// Compute a fine subdivision for the points.

		bKdTree.refine(splitRule);

		// Remove all objects but leave subdivision intact.

		bKdTree.eraseObjects();

		// Insert the first time-window.
		
		std::deque<ConstObjectIterator> bIteratorSet;
		for (integer i = 0;i < initialSamples;++i)
		{
			bIteratorSet.push_back(bKdTree.insert(&bPointSet[i]));
		}

		// Compute temporal mutual information
		// -----------------------------------

		integer tLeft = -sigma;
		integer tRight = sigma;

		std::deque<ConstObjectIterator>::iterator jointIter = 
			jointIteratorSet.begin();
		std::deque<ConstObjectIterator>::iterator aIter = 
			aIteratorSet.begin();
		std::deque<ConstObjectIterator>::iterator bIter = 
			bIteratorSet.begin();
		
		// For all time instants...
		for (integer t = 0;t < samples;++t)
		{
			searchAllNeighbors(
				jointKdTree,
				DepthFirst_SearchAlgorithm_PointKdTree(),
				randomAccessRange(jointIter, trials),
				kNearest - 1,
				kNearest,
				randomAccessRange(constantIterator(infinity<real>()), trials),
				0,
				normBijection,
				0,
				&distanceArray);

			real estimate = 0;

			std::vector<integer> countSet(trials, 0);

			countAllNeighbors(
				aKdTree,
				randomAccessRange(aIter, trials),
				randomAccessRange(distanceArray.begin(), trials),
				normBijection,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate -= digamma<real>(countSet[j]);
			}

			estimate /= trials;
			estimate += digamma<real>(kNearest);
			estimate += digamma<real>(samples);

			*result = estimate;
			++result;

			// Update time window
			std::advance(jointIter, trials);
			std::advance(aIter, trials);
			std::advance(bIter, trials);
			++tLeft;
			++tRight;

			if (tLeft > 0)
			{
				for (integer i = 0;i < trials;++i)
				{
					jointKdTree.erase(jointIteratorSet.front());
					jointIteratorSet.pop_front();
					
					aKdTree.erase(aIteratorSet.front());
					aIteratorSet.pop_front();
					
					bKdTree.erase(bIteratorSet.front());
					bIteratorSet.pop_front();
				}
			}
			if (tRight < samples)
			{
				integer index = tRight * trials;
				for (integer i = 0;i < trials;++i)
				{
					jointIteratorSet.push_back(
						jointKdTree.insert(&jointPointSet[index]));
					aIteratorSet.push_back(
						aKdTree.insert(&aPointSet[index]));
					bIteratorSet.push_back(
						bKdTree.insert(&bPointSet[index]));
					++index;
				}
			}
		}
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
