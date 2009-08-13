#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_KL_HPP

#include "tim/core/differential_entropy_kl.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/countingiterator.h>
#include <pastel/sys/randomaccessrange.h>

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

#include <algorithm>

namespace Tim
{

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		// This function computes the Kozachenko-Leonenko
		// estimator for the differential entropy.
		// Apparently the original paper seems to be:
		//
		// "A statistical estimate for the entropy of a random vector", 
		// Kozachenko, L.F. and Leonenko, N.N., (1987)
		// Problems Infor. Transmiss., 23(2), 9-16
		//
		// However, I couldn't find or read that paper online
		// so I cite:
		//
		// "Synchronization and Interdependence Measures
		// and their Applications to the Electroencephalogram
		// of Epilepsy Patients and Clustering of Data",
		// Alexander Kraskov, Ph.D. thesis, 2004

		const integer samples = signal->samples();
		const integer dimension = signal->dimension();

		Array<real, 2> distanceArray(1, samples);

		std::vector<PointD> pointSet;
		constructPointSet(signal, pointSet);

		typedef PointKdTree<real, Dynamic, Pointer_ObjectPolicy_PointKdTree<real, Dynamic> > KdTree;
		typedef typename KdTree::ConstObjectIterator ConstObjectIterator;
		
		KdTree kdTree(ofDimension(dimension));
		kdTree.insert(
			countingIterator(&pointSet.front()), 
			countingIterator(&pointSet.front() + pointSet.size()));

		std::vector<ConstObjectIterator> iteratorSet;
		iteratorSet.reserve(pointSet.size());
		std::copy(countingIterator(kdTree.begin()),
			countingIterator(kdTree.end()), 
			std::back_inserter(iteratorSet));

		kdTree.refine(SlidingMidpoint2_SplitRule_PointKdTree());

		searchAllNeighbors(
			kdTree,
			DepthFirst_SearchAlgorithm_PointKdTree(),
			randomAccessRange(iteratorSet.begin(), iteratorSet.end()),
			kNearest - 1,
			kNearest, 
			randomAccessRange(constantIterator(infinity<real>()), pointSet.size()),
			maxRelativeError,
			normBijection,
			0,
			&distanceArray);

		real estimate = 0;
#pragma omp parallel for reduction(+ : estimate)
		for (integer i = 0;i < samples;++i)
		{
			// Here we should add the logarithm of 
			// _twice_ the distance to the k:th neighbor.
			// However, we delay this by noting that:
			// log(dist * 2) = log(dist) + log(2)
			if (distanceArray(0, i) > 0)
			{
				estimate += normBijection.toLnNorm(distanceArray(0, i));
			}
		}
		// Here we take into account doubling the distances.
		estimate += samples * constantLn2<real>();

		estimate *= (real)dimension / samples;
		estimate -= digamma<real>(kNearest);
		estimate += digamma<real>(samples);
		// Here we add the logarithm of the volume of 
		// a sphere with _diameter_ 1. That is, the logarithm
		// of the volume of a sphere with radius 1/2:
		// log(unitVol * (1/2)^d) = log(unitVol) - d * log(2)
		estimate += normBijection.lnVolumeUnitSphere(dimension);
		estimate -= dimension * constantLn2<real>();

		return estimate;
	}

}

#endif
