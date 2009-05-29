#ifndef TIM_DIFFERENTIAL_ENTROPY_HPP
#define TIM_DIFFERENTIAL_ENTROPY_HPP

#include "tim/core/differential_entropy.h"

#include <pastel/geometry/all_nearest_neighbors_kdtree.h>

namespace Tim
{

	template <typename NormBijection>
	real differentialEntropy(
		const SignalPtr& signal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection)
	{
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

		const integer points = signal->size();
		const integer dimension = signal->dimension();

		Array<2, real> distanceArray(1, points);

		allNearestNeighborsKdTree(
			signal->pointSet(),
			kNearest - 1,
			kNearest,
			infinity<real>(),
			maxRelativeError,
			normBijection,
			0,
			&distanceArray);

		real estimate = 0;
		for (integer i = 0;i < points;++i)
		{
			// Here we should add the logarithm of 
			// _twice_ the distance to the k:th neighbor.
			// However, we delay this by noting that:
			// log(dist * 2) = log(dist) + log(2)
			estimate += normBijection.toLnNorm(distanceArray(0, i));
		}
		// Here we take into account doubling the distances.
		estimate += points * constantLn2<real>();

		estimate *= (real)dimension / points;
		estimate -= digamma<real>(kNearest);
		estimate += digamma<real>(points);
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
