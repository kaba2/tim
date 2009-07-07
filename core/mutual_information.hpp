#ifndef TIM_MUTUAL_INFORMATION_HPP
#define TIM_MUTUAL_INFORMATION_HPP

#include "tim/core/mutual_information.h"
#include "tim/core/signal_tools.h"

#include <pastel/geometry/search_all_neighbors_kdtree.h>
#include <pastel/geometry/count_all_neighbors_kdtree.h>
#include <pastel/geometry/distance_point_point.h>

namespace Tim
{

	template <typename NormBijection>
	real mutualInformation2(
		const SignalPtr& jointSignal,
		const std::vector<SignalPtr>& marginalSignalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		const integer signals = marginalSignalSet.size();
		const integer samples = jointSignal->samples();

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE2(jointSignal->samples() == marginalSignalSet[i]->samples(),
				jointSignal->samples(), marginalSignalSet[i]->samples());

			jointDimension += marginalSignalSet[i]->dimension();
		}

		ENSURE2(jointDimension == jointSignal->dimension(),
			jointDimension, jointSignal->dimension());

		Array<2, integer> nearestArray(1, samples);

		{
			std::vector<PointD> jointPointSet;
			constructPointSet(jointSignal, jointPointSet);

			searchAllNeighborsKdTree(
				jointPointSet,
				kNearest - 1,
				kNearest,
				infinity<real>(),
				maxRelativeError,
				normBijection,
				16,
				SlidingMidpoint2_SplitRule(),
				&nearestArray);
		}

		std::vector<real> distanceSet(samples);
		std::vector<integer> countSet(samples, 0);
		std::vector<PointD> marginalPointSet;

		real estimate = 0;
		real logVolumeSum = 0;

		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = marginalSignalSet[i];
			const integer marginalDimension = signal->dimension();
			
			constructPointSet(signal, marginalPointSet);

#pragma omp parallel for
			for (integer j = 0;j < samples;++j)
			{
				distanceSet[j] = 
					distance2(marginalPointSet[j],
					marginalPointSet[nearestArray(0, j)],
					normBijection);
			}

			countAllNeighborsKdTree(
				marginalPointSet,
				distanceSet,
				0,
				//maxRelativeError,
				normBijection,
				countSet);

			// We should actually add the logarithm of the
			// sphere with _diameter_ 1 (not radius 1). However, the radius
			// doesn't actually matter since it is cancelled
			// away later.
			logVolumeSum += 
				normBijection.lnVolumeUnitSphere(marginalDimension);

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= digamma<real>(countSet[j]);
			}
		}

		// Here the radii of the marginal spheres cancel
		// because the marginal dimensions sum to the joint dimension.
		// Note that with the infinity norm 'logVolumeSum' ends
		// up being zero. However, this is not the case with other
		// norms.
		logVolumeSum -= 
			normBijection.lnVolumeUnitSphere(jointSignal->dimension());

		estimate /= samples;
		estimate += logVolumeSum;
		estimate += digamma<real>(kNearest);
		estimate += (signals - 1) * digamma<real>(samples);
		estimate -= (real)(signals - 1) / kNearest;

		return estimate;
	}

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
		const integer samples = jointSignal->samples();

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE2(jointSignal->samples() == marginalSignalSet[i]->samples(),
				jointSignal->samples(), marginalSignalSet[i]->samples());

			jointDimension += marginalSignalSet[i]->dimension();
		}

		ENSURE2(jointDimension == jointSignal->dimension(),
			jointDimension, jointSignal->dimension());

		Array<2, real> distanceArray(1, samples);

		{
			std::vector<PointD> jointPointSet;
			constructPointSet(jointSignal, jointPointSet);

			searchAllNeighborsKdTree(
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

		real logVolumeSum = 0;

		for (integer i = 0;i < signals;++i)
		{
			const SignalPtr signal = marginalSignalSet[i];
			const integer marginalDimension = signal->dimension();
			
			std::vector<PointD> marginalPointSet;
			constructPointSet(signal, marginalPointSet);

			countAllNeighborsKdTree(
				marginalPointSet,
				distanceSet,
				0,
				//maxRelativeError,
				normBijection,
				countSet);

			// We should actually add the logarithm of the
			// sphere with _diameter_ 1 (not radius 1). However, the radius
			// doesn't actually matter since it is cancelled
			// away later.
			logVolumeSum += 
				normBijection.lnVolumeUnitSphere(marginalDimension);

#pragma omp parallel for reduction(+ : estimate)
			for (integer j = 0;j < samples;++j)
			{
				estimate -= digamma<real>(countSet[j]);
			}
		}

		// Here the radii of the marginal spheres cancel
		// because the marginal dimensions sum to the joint dimension.
		// Note that with the infinity norm 'logVolumeSum' ends
		// up being zero. However, this is not the case with other
		// norms.
		logVolumeSum -= 
			normBijection.lnVolumeUnitSphere(jointSignal->dimension());

		estimate /= samples;
		estimate += logVolumeSum;
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
		slice(jointSignal, marginalSignalSet);
		
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
		const SmallSet<integer>& partition,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		std::vector<SignalPtr> marginalSignalSet;
		slice(jointSignal, partition, marginalSignalSet);

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

	template <typename NormBijection>
	real mutualInformationFromEntropy(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
		ENSURE1(kNearest > 0, kNearest);
		ENSURE1(maxRelativeError >= 0, maxRelativeError);

		std::vector<SignalPtr> marginalSignalSet;
		slice(jointSignal, partition, marginalSignalSet);

		const integer signals = marginalSignalSet.size();
		const integer samples = jointSignal->samples();

		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			ENSURE2(jointSignal->samples() == marginalSignalSet[i]->samples(),
				jointSignal->samples(), marginalSignalSet[i]->samples());

			jointDimension += marginalSignalSet[i]->dimension();
		}

		ENSURE2(jointDimension == jointSignal->dimension(),
			jointDimension, jointSignal->dimension());

		real estimate = 0;

		for (integer i = 0;i < signals;++i)
		{
			estimate += differentialEntropy(
				marginalSignalSet[i],
				kNearest,
				maxRelativeError,
				normBijection);
		}

		estimate -= differentialEntropy(
			jointSignal,
			kNearest,
			maxRelativeError,
			normBijection);

		return estimate;
	}

}

#endif
