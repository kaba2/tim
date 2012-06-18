// TIM 1.2.0
// Kalle Rutanen
// http://kaba.hilvi.org
// Copyright (c) 2009 - 2011
//
// This library is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published 
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#ifndef TIM_EMBEDDING_FACTOR_HPP
#define TIM_EMBEDDING_FACTOR_HPP

#include "tim/core/embedding_factor_fn.h"

namespace Tim
{

	template <typename SignalPtr_Iterator>
	integer embeddingFactorFn(
		const ForwardIterator_Range<SignalPtr_Iterator>& signalSet)
	{
		typedef typename SignalPointSet::Point_ConstIterator
			Point_ConstIterator;

		if (signalSet.empty())
		{
			return 1;
		}

		// This function encapsulates the common
		// properties of the entropy estimation 
		// algorithms based on k-nearest neighbors.

		// This is done to avoid parallelization
		// issues with iterator range caching.

		signalSet.updateCache();

		SignalPointSet pointSet(
			signalSet, );

		const integer trials = signalSet.size();
		const integer samples = pointSet.samples();
		const integer dimension = signalSet.front()->dimension();
		const integer estimateSamples = samples * trials;

		Array<real> distanceArray(1, estimateSamples);

		// Find the distance to the k:th nearest neighbor for all points.

		searchAllNeighbors(
			pointSet.kdTree(),
			range(pointSet.begin(), pointSet.end()),
			kNearest - 1,
			kNearest, 
			(Array<Point_ConstIterator>*)0,
			&distanceArray,
			constantRange(infinity<real>(), estimateSamples),
			0,
			entropyAlgorithm.normBijection());

		// After we have found the distances, we simply evaluate
		// the generic entropy estimator over all samples.

		integer acceptedSamples = 0;
		real estimate = 0;
#pragma omp parallel for reduction(+ : estimate, acceptedSamples)
		for (integer i = 0;i < estimateSamples;++i)
		{
			// Points that are at identical positions do not
			// provide any information. Such samples are
			// not taken in the estimate.
			if (distanceArray(i) > 0)
			{
				estimate += entropyAlgorithm.sumTerm(distanceArray(i));
				++acceptedSamples;
			}
		}
		if (acceptedSamples > 0)
		{
			estimate = entropyAlgorithm.finishEstimate(
				estimate / acceptedSamples, dimension, kNearest, 
				estimateSamples);
		}
		else
		{
			// If all distances were zero, we can't say
			// anything about generic entropy. This is
			// marked with a NaN.
			estimate = nan<real>();
		}

		return estimate;
	}


}

#endif
