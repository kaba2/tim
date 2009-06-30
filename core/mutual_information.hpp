#ifndef TIM_MUTUAL_INFORMATION_HPP
#define TIM_MUTUAL_INFORMATION_HPP

#include "tim/core/mutual_information.h"

#include <pastel/geometry/all_nearest_neighbors_kdtree.h>
#include <pastel/geometry/count_all_nearest_neighbors_kdtree.h>
#include <pastel/geometry/count_all_nearest_neighbors_bruteforce.h>

namespace Tim
{

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& jointSignal,
		const std::vector<SignalPtr>& marginalSignalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		const integer signals = marginalSignalSet.size();
		const integer samples = jointSignal->height();

		/*
		real estimate2 = 0;

		for (integer i = 0;i < signals;++i)
		{
			estimate2 += differentialEntropy(
				marginalSignalSet[i],
				kNearest,
				maxRelativeError,
				normBijection);
		}

		estimate2 -= differentialEntropy(
			jointSignal,
			kNearest,
			maxRelativeError,
			normBijection);

		return estimate2;
		*/

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE2(jointSignal->height() == marginalSignalSet[i]->height(),
				jointSignal->height(), marginalSignalSet[i]->height());

			jointDimension += marginalSignalSet[i]->width();
		}

		ENSURE2(jointDimension == jointSignal->width(),
			jointDimension, jointSignal->width());

		Array<2, real> distanceArray(1, samples);

		{
			std::vector<PointD> jointPointSet;
			constructPointSet(jointSignal, jointPointSet);

			allNearestNeighborsKdTree(
				jointPointSet,
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

		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = marginalSignalSet[i];
			
			std::vector<PointD> marginalPointSet;
			constructPointSet(signal, marginalPointSet);

			countAllNearestNeighborsKdTree(
				marginalPointSet,
				distanceSet,
				0,
				//maxRelativeError,
				normBijection,
				countSet);

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= digamma<real>(countSet[j] + 1);
			}
		}

		estimate /= samples;
		estimate += digamma<real>(kNearest);
		estimate += (signals - 1) * digamma<real>(samples);

		return estimate;
	}

	template <typename NormBijection>
	real mutualInformation(
		const std::vector<SignalPtr>& marginalSignalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		const SignalPtr jointSignal =
			mergeSignalDimensions(marginalSignalSet);
		
		return Tim::mutualInformation(
			jointSignal,
			marginalSignalSet,
			kNearest,
			maxRelativeError,
			normBijection);
	}

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& jointSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		std::vector<SignalPtr> marginalSignalSet;
		splitDimensions(jointSignal, marginalSignalSet);
		
		return Tim::mutualInformation(
			jointSignal,
			marginalSignalSet,
			kNearest,
			maxRelativeError,
			normBijection);
	}

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& aMarginalSignal,
		const SignalPtr& bMarginalSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		std::vector<SignalPtr> marginalSignalSet;
		marginalSignalSet.push_back(aMarginalSignal);
		marginalSignalSet.push_back(bMarginalSignal);
		
		return Tim::mutualInformation(
			marginalSignalSet,
			kNearest,
			maxRelativeError,
			normBijection);
	}

}

#endif
