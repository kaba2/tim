#ifndef TIM_MUTUAL_INFORMATION_HPP
#define TIM_MUTUAL_INFORMATION_HPP

#include "tim/core/mutual_information.h"

#include <pastel/geometry/all_nearest_neighbors_kdtree.h>
#include <pastel/geometry/count_all_nearest_neighbors_kdtree.h>

namespace Tim
{

	template <typename NormBijection>
	real mutualInformation(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		const integer signals = signalSet.size();

		const SignalPtr jointSignal =
			mergeSignalDimensions(signalSet);

		const integer points = jointSignal->size();
		const integer jointDimension = jointSignal->dimension();

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
			const SignalPtr signal = signalSet[i];
			
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
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		std::vector<SignalPtr> signalSet;
		signalSet.push_back(aSignal);
		signalSet.push_back(bSignal);
		
		return Tim::mutualInformation(
			signalSet,
			kNearest,
			maxRelativeError,
			normBijection);
	}

}

#endif
