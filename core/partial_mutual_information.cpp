#include "tim/core/partial_mutual_information.h"

#include "tim/core/partial_mutual_information.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>
#include <pastel/geometry/distance_point_point.h>

#include <pastel/sys/math_functions.h>
#include <pastel/sys/pastelomp.h>

namespace Tim
{

	TIMCORE real partialMutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		const SignalPtr& cSignal,
		integer kNearest)
	{
		/*
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(aSignal->samples(), ==, bSignal->samples());
		ENSURE_OP(aSignal->samples(), ==, cSignal->samples());

		const real maxRelativeError = 0;
		const integer samples = aSignal->samples();

		// See
		// "Partial Mutual Information for Coupling Analysis 
		// of Multivariate Time Series", 
		// Stefan Frenzel and Bernd Pompe,
		// Physical Review Letters, 2007.

		const Infinity_NormBijection<real> normBijection;

		std::vector<SignalPtr> marginalSet;
		marginalSet.reserve(3);

		// We will define signals x, y, and z as
		// covariance-normalized signals a, b, and c.

		// We begin by forming a joint signal which
		// is formed by merging the samples from x, y, and z
		// together to form higher dimensional samples in xyz.

		// Note that the merging order of the signals
		// is important. This is because we need the
		// following subspaces:
		// xyz, xz, yz, and z. 
		// By merging the signals in the order (x, z, y),
		// we get the following subspaces by slicing:
		// xzy, xz, zy, and z.
		// Because merging order does not otherwise matter,
		// these will do fine.

		marginalSet.push_back(aSignal);
		marginalSet.push_back(cSignal);
		marginalSet.push_back(bSignal);
		
		const SignalPtr jointSignal = merge(
			marginalSet.begin(), marginalSet.end());
		
		const integer xBegin = 0;
		const integer xEnd = aSignal->dimension();
		const integer zBegin = xEnd;
		const integer zEnd = xEnd + cSignal->dimension();
		const integer yBegin = zEnd;
		const integer yEnd = zEnd + bSignal->dimension();

		const SignalPtr zSignal = slice(jointSignal, zBegin, zEnd);
		const SignalPtr xzSignal = slice(jointSignal, xBegin, zEnd);
		const SignalPtr yzSignal = slice(jointSignal, zBegin, yEnd);
				
		// For each sample point in the joint space,
		// find the distance to the k:th nearest neighbor.
		Array<2, real> distanceArray(1, samples);
		{
			std::vector<PointD> jointPointSet;
			constructPointSet(jointSignal, jointPointSet);

			searchAllNeighborsKdTree(
				jointPointSet,
				DepthFirst_SearchAlgorithm_PointKdTree(),
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				kNearest - 1,
				kNearest,
				infinity<real>(),
				maxRelativeError,
				normBijection,
				16,
				SlidingMidpoint2_SplitRule(),
				0,
				&distanceArray);
		}

		std::vector<real> distanceSet;
		distanceSet.reserve(samples);
		for (integer i = 0;i < samples;++i)
		{
			distanceSet.push_back(distanceArray(0, i));
		}

		real estimate = 0;

		std::vector<integer> countSet(samples, 0);

		// Orthogonally project the samples
		// into xz-subspace. For each sample, count the number of
		// neighbouring samples that fall in the
		// nn-distance just computed. Compute harmonic
		// number from that count and add to the estimate.
		{
			std::vector<PointD> marginalPointSet;
			constructPointSet(xzSignal, marginalPointSet);

			countAllNeighborsKdTree(
				marginalPointSet,
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				distanceSet.begin(),
				normBijection,
				16,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= harmonicNumber<real>(countSet[j]);
			}
		}

		// Do the same for the yz-subspace.
		{
			std::vector<PointD> marginalPointSet;
			constructPointSet(yzSignal, marginalPointSet);

			countAllNeighborsKdTree(
				marginalPointSet,
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				distanceSet.begin(),
				normBijection,
				16,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= harmonicNumber<real>(countSet[j]);
			}
		}

		// Do the same for the z-subspace.
		{
			std::vector<PointD> marginalPointSet;
			constructPointSet(zSignal, marginalPointSet);

			countAllNeighborsKdTree(
				marginalPointSet,
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				distanceSet.begin(),
				normBijection,
				16, 
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate += harmonicNumber<real>(countSet[j]);
			}
		}

		// Divide by the number of samples to get mean.
		estimate /= samples;

		// Finish off the estimate.
		estimate += harmonicNumber<real>(kNearest - 1);

		return estimate;
		*/
		return 0;
	}
		

}
