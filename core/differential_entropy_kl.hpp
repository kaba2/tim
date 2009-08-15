#ifndef TIM_DIFFERENTIAL_ENTROPY_KL_HPP
#define TIM_DIFFERENTIAL_ENTROPY_KL_HPP

#include "tim/core/differential_entropy_kl.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/countingiterator.h>
#include <pastel/sys/randomaccessrange.h>

#include <pastel/geometry/search_all_neighbors_pointkdtree.h>

#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>
#include <numeric>

namespace Tim
{

	template <
		typename Signal_Iterator, 
		typename Real_OutputIterator,
		typename NormBijection>
	void differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
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

		if (signalSet.empty())
		{
			return;
		}

		signalSet.updateCache();

		const integer trials = signalSet.size();
		const integer samples = minSamples(signalSet);
		const integer dimension = signalSet.front()->dimension();
		const integer timeWindowRadiusSamples = 2 * timeWindowRadius + 1;

		if (timeWindowRadius < samples - 1)
		{
			std::vector<real> estimateSet(samples);
#pragma omp parallel
			{
			Array<real, 2> distanceArray(1, trials);
			SignalPointSet pointSet(signalSet, 
				(timeWindowRadius < samples / 2) ? SignalPointSet_TimeWindow::StartEmpty : 
				SignalPointSet_TimeWindow::StartFull);

#pragma omp for
			for (integer t = 0;t < samples;++t)
			{
				const integer tLeft = std::max(t - timeWindowRadius, 0);
				const integer tRight = std::min(t + timeWindowRadius + 1, samples);
				const integer tDelta = t - tLeft;
				const integer tWidth = tRight - tLeft;

				pointSet.setTimeWindow(tLeft, tRight);
				
				searchAllNeighbors(
					pointSet.kdTree(),
					DepthFirst_SearchAlgorithm_PointKdTree(),
					randomAccessRange(pointSet.begin() + tDelta * trials, 
					pointSet.begin() + (tDelta + 1) * trials),
					kNearest - 1,
					kNearest, 
					randomAccessRange(constantIterator(infinity<real>()), trials),
					maxRelativeError,
					normBijection,
					0,
					&distanceArray);

				real estimate = 0;
				for (integer i = 0;i < trials;++i)
				{
					// Here we should add the logarithm of 
					// _twice_ the distance to the k:th neighbor.
					// However, we delay this by noting that:
					// log(dist * 2) = log(dist) + log(2)
					if (distanceArray(i) > 0)
					{
						estimate += normBijection.toLnNorm(distanceArray(i));
					}
				}
				// Here we take into account doubling the distances.
				estimate += trials * constantLn2<real>();

				estimate *= (real)dimension / trials;
				estimate -= digamma<real>(kNearest);
				estimate += digamma<real>(tWidth * trials);
				// Here we add the logarithm of the volume of 
				// a sphere with _diameter_ 1. That is, the logarithm
				// of the volume of a sphere with radius 1/2:
				// log(unitVol * (1/2)^d) = log(unitVol) - d * log(2)
				estimate += normBijection.lnVolumeUnitSphere(dimension);
				estimate -= dimension * constantLn2<real>();

				estimateSet[t] = estimate;
			}
			}

			std::copy(estimateSet.begin(), estimateSet.end(), result);
		}
		else
		{
			const real estimate = 
				Tim::differentialEntropy(
				signalSet, 
				kNearest, maxRelativeError, 
				normBijection);

			std::fill_n(result, samples, estimate);
		}
	}

	template <
		typename Signal_Iterator, 
		typename Real_OutputIterator>
	typename boost::disable_if<
		boost::is_convertible<Real_OutputIterator, integer> >::type
	differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest)
	{
		Tim::differentialEntropy(
			signalSet,
			timeWindowRadius,
			result,
			maxRelativeError,
			kNearest,
			Euclidean_NormBijection<real>());
	}

	template <
		typename Signal_Iterator, 
		typename NormBijection>
	real differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		const integer trials = signalSet.size();
		const integer samples = minSamples(signalSet);
		const integer dimension = signalSet.front()->dimension();
		const integer estimateSamples = samples * trials;

		Array<real, 2> distanceArray(1, estimateSamples);
		SignalPointSet pointSet(signalSet, SignalPointSet_TimeWindow::StartFull);
	
		searchAllNeighbors(
			pointSet.kdTree(),
			DepthFirst_SearchAlgorithm_PointKdTree(),
			randomAccessRange(pointSet.begin(), pointSet.end()),
			kNearest - 1,
			kNearest, 
			randomAccessRange(constantIterator(infinity<real>()), estimateSamples),
			maxRelativeError,
			normBijection,
			0,
			&distanceArray);

		real estimate = 0;
#pragma omp parallel for reduction(+ : estimate)
		for (integer i = 0;i < estimateSamples;++i)
		{
			// Here we should add the logarithm of 
			// _twice_ the distance to the k:th neighbor.
			// However, we delay this by noting that:
			// log(dist * 2) = log(dist) + log(2)
			if (distanceArray(i) > 0)
			{
				estimate += normBijection.toLnNorm(distanceArray(i));
			}
		}
		// Here we take into account doubling the distances.
		estimate += estimateSamples * constantLn2<real>();

		estimate *= (real)dimension / estimateSamples;
		estimate -= digamma<real>(kNearest);
		estimate += digamma<real>(estimateSamples);
		// Here we add the logarithm of the volume of 
		// a sphere with _diameter_ 1. That is, the logarithm
		// of the volume of a sphere with radius 1/2:
		// log(unitVol * (1/2)^d) = log(unitVol) - d * log(2)
		estimate += normBijection.lnVolumeUnitSphere(dimension);
		estimate -= dimension * constantLn2<real>();

		return estimate;
	}

	template <typename Signal_Iterator>
	real differentialEntropy(
		const ForwardRange<Signal_Iterator>& signalSet,
		real maxRelativeError,
		integer kNearest)
	{
		return Tim::differentialEntropy(
			signalSet,
			maxRelativeError,
			kNearest,
			Euclidean_NormBijection<real>());
	}

	// Overloads for a single signal
	// -----------------------------

	template <
		typename Real_OutputIterator,
		typename NormBijection>
	void differentialEntropy(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(timeWindowRadius, >=, 0);
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		Tim::differentialEntropy(
			forwardRange(constantIterator(signal)),
			timeWindowRadius,
			result,
			maxRelativeError,
			kNearest,
			normBijection);
	}

	template <typename Real_OutputIterator>
	typename boost::disable_if<
		boost::is_convertible<Real_OutputIterator, integer> >::type
	differentialEntropy(
		const SignalPtr& signal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		real maxRelativeError,
		integer kNearest)
	{
		Tim::differentialEntropy(
			signal,
			timeWindowRadius,
			result,
			maxRelativeError,
			kNearest, 
			Euclidean_NormBijection<real>());
	}

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		real maxRelativeError,
		integer kNearest,
		const NormBijection& normBijection)
	{
		ENSURE_OP(kNearest, >, 0);
		ENSURE_OP(maxRelativeError, >=, 0);

		return Tim::differentialEntropy(
			forwardRange(constantIterator(signal)),
			maxRelativeError,
			kNearest,
			normBijection);
	}

}

#endif
