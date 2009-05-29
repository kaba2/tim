#ifndef TIM_MUTUAL_INFORMATION_H
#define TIM_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

namespace Tim
{

	//! Computes mutual information between a set of signals.
	/*!
	Preconditions:
	jointSignal->size() == marginalSignalSet[i]->size()
	jointSignal->dimension() = sum_i marginalSignalSet[i]->dimension()
	kNearest > 0
	maxRelativeError >= 0
	
	This is the most general version of the functions
	to compute mutual information. Here by 'mutual information'
	we actually mean 'total correlation' which is one generalization
	of mutual information for higher number of signals than two.
	If H(X) denotes the differential entropy of a real 
	random variable X, and we are given n random variables
	(X_1, ..., X_n)	then the total correlation between 
	them is defined by:
	I(X_1, ..., X_n) = (sum_i H(X_i)) - H(X_1, ..., X_n)
	While differential entropy does not make
	sense w.r.t the amount of information, total correlation does.
	See 'tim/core/differential_entropy.h' for more information
	on differential entropy.
	
	The mutual information is estimated by an algorithm
	that uses nearest neighbor searching heavily.

	jointSignal:
	The joint signal must contain the same data 
	as the marginal signals, arranged subsequently
	in subdimension ranges. For example, 
	dimensions 0 contain 1-d marginalSignalSet[0],
	dimensions 1-2 contains 2-d marginalSignalSet[1], and
	dimension 3 contains 1-d marginalSignalSet[2].
	
	marginalSignalSet:
	A set of marginal signals. They must all
	have the same size. However, their dimension
	can vary.

	normBijection:
	Different norms can be used in the nearest
	neighbor search. The norms are varied by supplying
	different kinds of 'norm bijections'. 
	See 'pastel/math/normbijection.h' in the Pastel library
	for documentation and predefined norm bijections.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	maxRelativeError:
	The maximum relative error in the distance of the
	nearest neighbors. Allowing approximate results in
	nearest neighbor searching can accelerate performance 
	by one to two orders of magnitude in higher dimensions.
	0 = exact, 1 = 100% relative error, 2 = 200% relative error,
	etc.
	*/

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& jointSignal,
		const std::vector<SignalPtr>& marginalSignalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	//! Computes mutual information between marginal signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation() by first
	computing the joint signal from the marginal signals.
	See the documentation for the more general mutualInformation().
	*/

	template <typename NormBijection>
	real mutualInformation(
		const std::vector<SignalPtr>& marginalSignalSet,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	//! Computes mutual information between 1-d marginal signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation() by splitting
	1-dimensional marginal signals from the joint signal.
	Splitting does not require copying: the
	same data is aliased.
	See the documentation for the more general mutualInformation().
	*/

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& jointSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	//! Computes mutual information between two marginal signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation() by pushing
	the marginal signals into a vector.
	See the documentation for the more general mutualInformation().
	*/

	template <typename NormBijection>
	real mutualInformation(
		const SignalPtr& aMarginalSignal,
		const SignalPtr& bMarginalSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection);

	TIMCORE real correlatedGaussianMutualInformation(
		const DynamicMatrix& correlation);

}

#include "tim/core/mutual_information.hpp"

#endif
