#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/math/normbijection.h>

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>

namespace Tim
{

	TIMCORE void transferEntropy(
		const std::vector<SignalPtr>& aEnsemble,
		const std::vector<SignalPtr>& aFutureEnsemble,
		const std::vector<SignalPtr>& bEnsemble,
		const Array<2, SignalPtr>& cEnsembleSet,
		integer sigma,
		integer kNearest,
		std::vector<real>& estimateSet)
	{
		ENSURE_OP(sigma, >=, 0);
		ENSURE_OP(kNearest, >, 0);

		const integer trials = aEnsemble.size();

		// Check that each ensemble contains the same number
		// of trials.

		ENSURE_OP(aFutureEnsemble.size(), ==, trials);
		ENSURE_OP(bEnsemble.size(), ==, trials);
		ENSURE_OP(cEnsembleSet.width(), ==, trials);

		if (trials == 0)
		{
			return;
		}

		const integer cSize = cEnsembleSet.height();

		// 1) Check that ensemble signals are of the same
		// dimension. 
		// 2) Because the delay-embeddings can result in signals
		// with slightly different number of samples,
		// we shall use the minimum of the sample count
		// among the signals. So we shall also find
		// the minimum sample count.

		const integer aDimension = aEnsemble.front()->dimension();
		const integer aFutureDimension = aFutureEnsemble.front()->dimension();
		const integer bDimension = bEnsemble.front()->dimension();
		const integer cDimension = cSize > 0 ? cEnsembleSet(0, 0)->dimension() : 0;

		integer samples = 0;

		for (integer i = 0;i < trials;++i)
		{
			ENSURE_OP(aDimension, ==, aEnsemble[i]->dimension());
			ENSURE_OP(bDimension, ==, bEnsemble[i]->dimension());
			ENSURE_OP(aFutureDimension, ==, aFutureEnsemble[i]->dimension());

			samples = std::min(samples, aEnsemble[i]->samples());
			samples = std::min(samples, aFutureEnsemble[i]->samples());
			samples = std::min(samples, bEnsemble[i]->samples());

			for (integer j = 0;j < cSize;++j)
			{
				ENSURE_OP(cDimension, ==, cEnsembleSet(i, j)->dimension());

				samples = std::min(samples, cEnsembleSet(i, j)->samples());
			}
		}

		// Find the sum of the dimensions of the cEnsemble.

		integer cDimensionTotal = 0;
		for (integer j = 0;j < cSize;++j)
		{
			cDimensionTotal += cEnsembleSet(0, j)->dimension();
		}

		// Form joint signals for each trial.

		// We need the following joint signals
		// (w is the future of X):
		//
		// 1) wXYZ
		// 2) XYZ
		// 3) XZ
		// 4) wXz
		//
		// By a clever reordering, we get away by just forming
		// the biggest joint signal wXZY and aliasing the
		// other joint signals from that:
		//
		// 1) wXZY
		// 2)  XZY
		// 3)  XZ
		// 4) wXZ

		std::vector<SignalPtr> jointEnsemble;
		jointEnsemble.reserve(trials);

		for (integer i = 0;i < trials;++i)
		{
			std::vector<SignalPtr> jointSet;
			jointSet.reserve(3 + cSize);

			jointSet.push_back(aEnsemble[i]);
			jointSet.push_back(aFutureEnsemble[i]);
			for (integer j = 0;j < cSize;++j)
			{
				jointSet.push_back(cEnsembleSet(i, j));
			}
			jointSet.push_back(bEnsemble[i]);

			jointEnsemble.push_back(merge(jointSet));
		}

		/*
		From now on, to understand what we are doing, 
		see the papers:

		"Assessing Coupling Dynamics from an Ensemble of Multivariate
		Time-series.", 
		German Gómez Herrero et al.
		
		"Measuring Information Transfer", 
		Thomas Schreiber,
		Physical Review Letters,
		vol. 85, num. 2, 2000.
		*/

		// Slice the smaller joint signals from the
		// bigger one. Note that here we need the signals 
		// to be merged in the right order.

		const integer wBegin = 0;
		const integer wEnd = aFutureDimension;
		const integer xBegin = wEnd;
		const integer xEnd = xBegin + aDimension;
		const integer zBegin = xEnd;
		const integer zEnd = xBegin + cDimensionTotal;
		const integer yBegin = zEnd;
		const integer yEnd = yBegin + bDimension;

		const integer xzyDimension = yEnd - xBegin;
		const integer xzDimension = zEnd - xBegin;
		const integer wxzDimension = zEnd - wBegin;
		const integer jointDimension = yEnd - wBegin;

		const Infinity_NormBijection<real> normBijection;
		Array<2, real> distanceArray(1, samples);

		const integer sigmaSamples = 2 * sigma + 1;

		std::vector<integer> countSet(trials, 0);
		std::vector<real> distanceSet(samples, 0);
		std::vector<PointD> pointSet;
		pointSet.reserve(trials * sigmaSamples);

		estimateSet.resize(samples);
		std::fill(estimateSet.begin(), estimateSet.end(), 0);

		for (integer i = sigma;i < samples - sigma;++i)
		{
			const integer sampleBegin = i - sigma;
			const integer sampleEnd = i + sigma + 1;

			// Compute for the current point its distance 
			// to the k:th nearest neighbor in the joint space,
			// which includes all samples from the ensemble
			// across the time window of sigma width.

			constructPointSet(jointEnsemble, 
				sampleBegin, sampleEnd, 
				wBegin, yEnd, pointSet);

			PointD jointPoint(
				ofDimension(jointDimension),
				withAliasing((real*)0));

			const ConstSparseIterator<CountingIterator<integer> >
				sparseIndexBegin(CountingIterator<integer>(sigma), sigmaSamples);

			searchAllNeighborsKdTree(
				pointSet,
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

			std::copy(
				distanceArray.columnBegin(0),
				distanceArray.columnEnd(0),
				distanceSet.begin());

			// Project all points to wXZ and count the number of
			// points that are within the computed nn-distance.

			constructPointSet(jointEnsemble,
				sampleBegin, sampleEnd, 
				wBegin, zEnd, pointSet);

			countAllNeighborsKdTree(
				pointSet,
				sparseIndexBegin,
				sparseIndexBegin + trials,
				distanceSet.begin(),
				normBijection,
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
				distanceSet.begin(),
				normBijection,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate += harmonicNumber<real>(countSet[j]);
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
				distanceSet.begin(),
				normBijection,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < trials;++j)
			{
				estimate -= harmonicNumber<real>(countSet[j]);
			}

			estimate /= trials;
			estimate -= harmonicNumber<real>(kNearest - 1);

			estimateSet[i] = estimate;
		}
	}

}
