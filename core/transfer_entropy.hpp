#ifndef TIM_TRANSFER_ENTROPY_HPP
#define TIM_TRANSFER_ENTROPY_HPP

#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/math/normbijection.h>

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>

namespace Tim
{

	template <typename Signal_ForwardIterator>
	void transferEntropy(
		const Array<2, SignalPtr>& embeddedSet,
		integer xIndex,
		integer yIndex,
		const Signal_ForwardIterator& xFutureBegin,
		const Signal_ForwardIterator& xFutureEnd,
		integer sigma,
		integer kNearest,
		std::vector<real>& estimateSet)
	{
		// See the paper:
		// "Assessing Coupling Dynamics from an Ensemble of Multivariate
		// Time-series.", 
		// German Gómez Herrero et al.

		const integer trials = embeddedSet.width();
		const integer signals = embeddedSet.height();

		ENSURE_OP(xIndex, >=, 0);
		ENSURE_OP(xIndex, <, signals);
		ENSURE_OP(yIndex, >=, 0);
		ENSURE_OP(yIndex, <, signals);
		ENSURE_OP(sigma, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(std::distance(xFutureBegin, xFutureEnd), ==, trials);
		ENSURE(equalDimension(xFutureBegin, xFutureEnd));

		const integer samples = 
			std::min(minSamples(xFutureBegin, xFutureEnd),
			minSamples(embeddedSet.begin(), embeddedSet.end()));

		if (trials == 0 || signals == 0 || samples == 0)
		{
			return;
		}
	
		// Form joint signals for each trial.

		// We need the following joint signals
		// (w is the future of x):
		//
		// 1) wXYZ
		// 2) XYZ
		// 3) XZ
		// 4) wXz
		//
		// By reordering, we can just form the biggest
		// joint signal wXZY and use its memory for
		// all point sets:
		//
		// 1) wXZY
		// 2)  XZY
		// 3)  XZ
		// 4) wXZ

		std::vector<SignalPtr> jointEnsemble;
		jointEnsemble.reserve(trials);

		Signal_ForwardIterator xFutureIter = xFutureBegin;
		for (integer i = 0;i < trials;++i)
		{
			std::vector<SignalPtr> jointSet;
			jointSet.reserve(1 + signals);

			jointSet.push_back(*xFutureIter);
			jointSet.push_back(embeddedSet(i, xIndex));

			for (integer j = 0;j < signals;++j)
			{
				if (j != xIndex && j != yIndex)
				{
					jointSet.push_back(embeddedSet(i, j));
				}
			}
			jointSet.push_back(embeddedSet(i, yIndex));

			jointEnsemble.push_back(
				merge(jointSet.begin(), jointSet.end()));

			++xFutureIter;
		}

		// Compute the dimensions of the signals.

		const integer jointDimension = jointEnsemble.front()->dimension();
		const integer xDimension = embeddedSet(0, xIndex)->dimension();
		const integer yDimension = embeddedSet(0, yIndex)->dimension();
		const integer zDimension = jointDimension - (xDimension * 2 + yDimension);
		const integer xFutureDimension = xDimension;

		const integer wBegin = 0;
		const integer wEnd = xFutureDimension;
		const integer xBegin = wEnd;
		const integer xEnd = xBegin + xDimension;
		const integer zBegin = xEnd;
		const integer zEnd = zBegin + zDimension;
		const integer yBegin = zEnd;
		const integer yEnd = yBegin + yDimension;

		const integer sigmaSamples = 2 * sigma + 1;

		// Data structures for nearest neighbors searching.

		const Infinity_NormBijection<real> normBijection;
		Array<2, real> distanceArray(1, trials);
		std::vector<PointD> pointSet;
		pointSet.reserve(trials * sigmaSamples);

		// Data structures for nearest neighbors counting.

		std::vector<integer> countSet(trials, 0);

		estimateSet.resize(samples);
		std::fill(estimateSet.begin(), estimateSet.end(), 0);

		integer percent = 0;
		for (integer i = sigma;i < samples - sigma;++i)
		{
			integer newPercent = (real)(i * 100) / samples;
			if (newPercent != percent)
			{
				log() << newPercent << "%, ";
				percent = newPercent;
			}

			const integer sampleBegin = i - sigma;
			const integer sampleEnd = i + sigma + 1;

			// Compute for the current point its distance 
			// to the k:th nearest neighbor in the joint space,
			// which includes all samples from the ensemble
			// across the time window of sigma radius.

			constructPointSet(jointEnsemble, 
				sampleBegin, sampleEnd, 
				wBegin, yEnd, pointSet);

			const ConstSparseIterator<CountingIterator<integer> >
				sparseIndexBegin(CountingIterator<integer>(sigma), sigmaSamples);

			searchAllNeighborsKdTree(
				pointSet,
				DepthFirst_SearchAlgorithm_PointKdTree(),
				sparseIndexBegin,
				sparseIndexBegin + trials,
				kNearest - 1,
				kNearest,
				infinity<real>(),
				0,
				normBijection,
				16,
				SlidingMidpoint2_SplitRule(),
				0,
				&distanceArray);

			// Project all points to wXZ and count the number of
			// points that are within the computed nn-distance.

			constructPointSet(jointEnsemble,
				sampleBegin, sampleEnd, 
				wBegin, zEnd, pointSet);

			countAllNeighborsKdTree(
				pointSet,
				sparseIndexBegin,
				sparseIndexBegin + trials,
				distanceArray.begin(),
				normBijection,
				16, 
				countSet.begin());

			real estimate = 0;

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate -= harmonicNumber<real>(countSet[j]);
			}

			// Project all points to XZY and count the number of
			// points that are within the computed nn-distance.

			constructPointSet(jointEnsemble, 
				sampleBegin, sampleEnd, 
				xBegin, yEnd, pointSet);

			countAllNeighborsKdTree(
				pointSet,
				sparseIndexBegin,
				sparseIndexBegin + trials,
				distanceArray.begin(),
				normBijection,
				16,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate -= harmonicNumber<real>(countSet[j]);
			}

			// Project all points to XZ and count the number of
			// points that are within the computed nn-distance.

			constructPointSet(jointEnsemble, 
				sampleBegin, sampleEnd, 
				xBegin, zEnd, pointSet);

			countAllNeighborsKdTree(
				pointSet,
				sparseIndexBegin,
				sparseIndexBegin + trials,
				distanceArray.begin(),
				normBijection,
				16,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate += harmonicNumber<real>(countSet[j]);
			}

			estimate /= trials;
			estimate += harmonicNumber<real>(kNearest - 1);

			estimateSet[i] = estimate;
		}
	}

}

#endif
