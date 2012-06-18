// Description: Local estimator concept
// Documentation: localestimators.txt

#ifndef TIM_LOCALESTIMATOR_CONCEPT_H
#define TIM_LOCALESTIMATOR_CONCEPT_H

#include "tim/core/mytypes.h"

namespace Tim
{

	namespace LocalEstimator_Concept
	{

		class LocalEstimator
		{
		public:
			class Instance
			{
			public:
				//! Constructs an instance of the local estimator.
				/*
				k:
				The k to use to search for the k:th nearest neighbor.
				n:
				The number of points in the sample.
				*/
				Instance(
					integer k,
					integer n);

				//! Computes the local joint estimate.
				/*!
				The local joint estimate is the part of the
				estimate contributed by the probability mass 
				in the k-nearest neighbor ball.
				*/
				real localJointEstimate() const;

				//! Computes the local marginal estimate.
				/*!
				The local marginal estimate is the part of the
				estimate contributed by the probability mass 
				in a marginal projection of the k-nearest neighbor 
				ball. When the points are considered by their
				marginal projections, the marginal projection of
				the k-nearest neighbor ball contains kMarginal
				number of points.
				*/
				real localMarginalEstimate(integer kMarginal) const;
			};
		};

	}

}

#endif
