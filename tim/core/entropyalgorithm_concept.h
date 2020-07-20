// Description: EntropyAlgorithm concept
// Documentation: generic_entropy.txt

#ifndef TIM_ENTROPYALGORITHM_CONCEPT_H
#define TIM_ENTROPYALGORITHM_CONCEPT_H

#include "tim/core/mytypes.h"

#include "pastel/math/normbijection/normbijection_concept.h"

namespace Tim
{

	class EntropyAlgorithm
	{
	public:
		// Defines the norm bijection to use
		// for k-nearest neighbor searching.
		typedef UserDefinedType Norm;
		
		// The norm bijection object is stored in 
		// the entropy algorithm object so that
		// it can be initialized properly
		// (e.g. Minkowski_Norm).
		const Norm& norm() const;
		
		// Compute a term in the sum based
		// on the distance to another point.
		// Note: the distance is in terms of the
		// norm.
		dreal sumTerm(Distance_Concept auto distance) const;

		// Apply some final transformation to
		// the estimate before being stored.
		dreal finishEstimate(
			dreal estimate, 
			integer dimension, 
			integer kNearest, 
			integer estimateSamples) const;
	};

}

#endif
