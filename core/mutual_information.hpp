#ifndef TIM_MUTUAL_INFORMATION_HPP
#define TIM_MUTUAL_INFORMATION_HPP

#include "tim/core/mutual_information.h"

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
			std::vector<DynamicPoint> pointSet;
			constructPointSet(jointSignal->view(),
				pointSet);

			allNearestNeighborsKdTree(
				pointSet,
				kNearest - 1,
				kNearest,
				infinity<real>(),
				maxRelativeError,
				normBijection,
				0,
				&distanceArray);
		}

		/*
		for (integer i = 0;i < points;++i)
		{
			distanceArray(0, i) = normBijection.toNorm(
				distanceArray(0, i));
		}
		*/
		
		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = signalSet[i];
			
			std::vector<DynamicPoint> pointSet;
			constructPointSet(signal->view(), pointSet);

			countAllNearestNeighborsKdTree(
				pointSet,
				0,
				pointSet.size() - 1,
				distanceArray,
				maxRelativeError,
				normBijection);
		}
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
