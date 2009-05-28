#ifndef TIM_MUTUAL_INFORMATION_HPP
#define TIM_MUTUAL_INFORMATION_HPP

#include "tim/core/mutual_information.h"

#include <pastel/geometry/all_nearest_neighbors_kdtree.h>
#include <pastel/geometry/count_all_nearest_neighbors_kdtree.h>

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
		const integer points = jointSignal->size();
		
		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE1(jointSignal->size() == marginalSignalSet[i]->size()
				jointSignal->size(), marginalSignalSet[i]->size());

			jointDimension += marginalSignalSet[i]->dimension();
		}

		ENSURE2(jointDimension == jointSignal->dimension(),
			jointDimension, jointSignal->dimension());

		Array<2, real> distanceArray(1, points);

		{
			allNearestNeighborsKdTree(
				jointSignal->pointSet(),
				kNearest - 1,
				kNearest,
				infinity<real>(),
				maxRelativeError,
				normBijection,
				0,
				&distanceArray);
		}

		real estimate = 0;
		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = marginalSignalSet[i];
			
			const integer totalNeighbors =
				countAllNearestNeighborsKdTree(
				signal->pointSet(),
				0,
				signal->size() - 1,
				distanceArray,
				maxRelativeError,
				normBijection);

			// Should this be digamma<real>(totalNeighbors) or 
			// digamma<real>(totalNeighbors + 1)?
			estimate -= digamma<real>(totalNeighbors + 1);
		}

		estimate /= points;
		estimate += digamma<real>(kNearest);
		estimate += (dimension - 1) * digamma<real>(points);

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
		
		return Pastel::mutualInformation(
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
		
		return Pastel::mutualInformation(
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
