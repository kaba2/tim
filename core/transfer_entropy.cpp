#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/math/normbijection.h>

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>

namespace Tim
{

	TIMCORE real transferEntropy(
		const SignalPtr& aEmbedded,
		const SignalPtr& aFuture,
		const SignalPtr& bEmbedded,
		const std::vector<SignalPtr>& cEmbeddedSet,
		integer sigma,
		integer kNearest)
	{
		ENSURE((aEmbedded->dimension() % aFuture->dimension()) == 0);
		ENSURE(!aEmbedded.empty());
		ENSURE(!aFuture.empty());
		ENSURE(!bEmbedded.empty());

		// Because the delay-embeddings can result in signals
		// with slightly different number of samples,
		// we shall use the minimum of the sample count
		// among the signals.

		// Find the minimum sample count.

		integer samples =
			std::min(aEmbedded->samples(),
			std::min(aFuture->samples(), 
			bEmbedded->samples()));

		integer cDimension = 0;
		const integer cSize = cEmbeddedSet.size();
		for (integer i = 0;i < cSize;++i)
		{
			const SignalPtr cEmbedded = cEmbeddedSet[i];

			ENSURE(!cEmbedded.empty());

			const integer cSamples = cEmbedded->samples();
			if (cSamples < samples)
			{
				samples = cSamples;
			}
			
			cDimension += cEmbedded->dimension();
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

		// We need the following joint signals
		// (w is the future of X):
		//
		// 1) wXYZ
		// 2) XYZ
		// 3) XZ
		// 4) wXz
		//
		// By a clever reordering, we get away by just forming
		// the biggest joint signal wXZY and slicing the
		// other joint signals from that:
		//
		// 1) wXZY
		// 2)  XZY
		// 3)  XZ
		// 4) wXZ

		// Form the joint signal.

		std::vector<SignalPtr> jointSet;
		jointSet.reserve(3 + cSize);
		jointSet.push_back(aFuture);
		jointSet.push_back(aEmbedded);
		for (integer i = 0;i < cSize;++i)
		{
			jointSet.push_back(cEmbeddedSet[i]);
		}
		jointSet.push_back(bEmbedded);

		const SignalPtr jointSignal = merge(jointSet);

		// Slice the smaller joint signals from
		// bigger one. Note here we need the signals to be merged
		// in the right order.

		const integer wBegin = 0;
		const integer wEnd = aFuture->dimension();
		const integer xBegin = wEnd;
		const integer xEnd = xBegin + aEmbedded->dimension();
		const integer zBegin = xEnd;
		const integer zEnd = xBegin + cDimension;
		const integer yBegin = zEnd;
		const integer yEnd = yBegin + bEmbedded->dimension();

		const SignalPtr xzySignal = slice(jointSignal, xBegin, yEnd);
		const SignalPtr xzSignal = slice(jointSignal, xBegin, zEnd);
		const SignalPtr wxzSignal = slice(jointSignal, wBegin, zEnd);

		const integer xzyDimension = xzySignal->dimension();
		const integer xzDimension = xzSignal->dimension();
		const integer wxzDimension = wxzSignal->dimension();

		// Compute for each point in the joint space its
		// distance to the k:th nearest neighbor.

		const InfinityNormBijection<real> normBijection;
		Array<2, real> distanceArray(1, samples);

		std::vector<PointD> jointPointSet;
		constructPointSet(jointSignal, jointPointSet);

		searchAllNeighborsKdTree(
			jointPointSet,
			kNearest - 1,
			kNearest,
			infinity<real>(),
			0,
			normBijection,
			16,
			SlidingMidpoint2_SplitRule(),
			0,
			&distanceArray);

		std::vector<real> distanceSet;
		distanceSet.reserve(samples);
		for (integer i = 0;i < samples;++i)
		{
			distanceSet.push_back(distanceArray(0, i));
		}

		std::vector<integer> countSet(samples);
		std::vector<real> estimate(samples, 0);

		const integer sigmaSamples = 2 * sigma + 1;

		std::vector<PointD> marginalPointSet;

		// Project all points to wXZ and count the number of
		// points that are within the computed nn-distance.

		for (integer i = sigma;i < samples - sigma;++i)
		{
			constructPointSet(jointSignal, wBegin, zEnd, marginalPointSet);

			countAllNeighborsKdTree(
				marginalPointSet,
				distanceSet,
				0,
				normBijection,
				countSet);

#pragma omp parallel for
			for (integer j = 0;j < samples;++j)
			{
				estimate[j] -= digamma<real>(countSet[j]);
			}
						
		}

	}

}
