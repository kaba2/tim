#include "tim/core/mutual_information.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/string_tools.h>

#include <pastel/math/matrix_tools.h>
#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/normbijection.h>

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>

#include <pastel/gfx/pcx.h>
#include <pastel/gfx/image_tools.h>

namespace Tim
{

	TIMCORE real mutualInformation(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		real maxRelativeError)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		if (signalSet.empty())
		{
			return 0;
		}

		const integer signals = signalSet.size();

		const SignalPtr jointSignal = merge(signalSet);
		const integer jointDimension = jointSignal->dimension();
		const integer samples = jointSignal->samples();

		std::vector<PointD> pointSet;
		constructPointSet(jointSignal, pointSet);

		const Infinity_NormBijection<real> normBijection;
		Array<2, real> distanceArray(1, samples);

		searchAllNeighborsKdTree(
			pointSet,
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

		std::vector<real> distanceSet;
		distanceSet.reserve(samples);

		std::copy(
			distanceArray.columnBegin(0),
			distanceArray.columnEnd(0),
			std::back_inserter(distanceSet));

		real estimate = 0;

		std::vector<integer> countSet(samples, 0);
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			const integer dimension = signalSet[i]->dimension();
			
			constructPointSet(
				jointSignal,
				0, samples,
				dimensionOffset, 
				dimensionOffset + dimension,
				pointSet);

			countAllNeighborsKdTree(
				pointSet,
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				distanceSet.begin(),
				normBijection,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= digamma<real>(countSet[j]);
			}

			dimensionOffset += dimension;
		}

		estimate /= samples;
		estimate += digamma<real>(kNearest);
		estimate += (signals - 1) * digamma<real>(samples);

		return estimate;
	}

	TIMCORE real mutualInformation2(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		real maxRelativeError)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		const integer signals = signalSet.size();
		const integer samples = signalSet.front()->samples();

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE_OP(signalSet[i]->samples(), ==, samples);

			jointDimension += signalSet[i]->dimension();
		}

		const Infinity_NormBijection<real> normBijection;

		const SignalPtr jointSignal = merge(signalSet);

		Array<2, integer> nearestArray(1, samples);

		std::vector<PointD> pointSet;
		constructPointSet(jointSignal, pointSet);

		searchAllNeighborsKdTree(
			pointSet,
			CountingIterator<integer>(0),
			CountingIterator<integer>(samples),
			kNearest - 1,
			kNearest,
			infinity<real>(),
			maxRelativeError,
			normBijection,
			16,
			SlidingMidpoint2_SplitRule(),
			&nearestArray);

		std::vector<real> distanceSet(samples);
		std::vector<integer> countSet(samples, 0);

		real estimate = 0;
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			const integer dimension = signalSet[i]->dimension();
			
			constructPointSet(
				jointSignal,
				0, samples,
				dimensionOffset,
				dimensionOffset + dimension,
				pointSet);

#pragma omp parallel for
			for (integer j = 0;j < samples;++j)
			{
				distanceSet[j] = 
					distance2(pointSet[j],
					pointSet[nearestArray(0, j)],
					normBijection);
			}

			countAllNeighborsKdTree(
				pointSet,
				CountingIterator<integer>(0),
				CountingIterator<integer>(samples),
				distanceSet.begin(),
				normBijection,
				countSet.begin());

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= digamma<real>(countSet[j]);
			}

			dimensionOffset += dimension;
		}

		estimate /= samples;
		estimate += digamma<real>(kNearest);
		estimate += (signals - 1) * digamma<real>(samples);
		estimate -= (real)(signals - 1) / kNearest;

		return estimate;
	}

	TIMCORE real mutualInformation(
		const SignalPtr& jointSignal,
		integer kNearest,
		real maxRelativeError)
	{
		std::vector<SignalPtr> signalSet;
		slice(jointSignal, signalSet);
		
		return Tim::mutualInformation(
			signalSet,
			kNearest,
			maxRelativeError);
	}

}

