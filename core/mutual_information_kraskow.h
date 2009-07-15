// Description: Mutual information estimation via Kraskow's algorithm

#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_H
#define TIM_MUTUAL_INFORMATION_KRASKOW_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/smallset.h>

namespace Tim
{

	//! Computes mutual information between a set of signals.
	/*!
	Preconditions:
	jointSignal->samples() == signalSet[i]->samples()
	jointSignal->dimension() = sum_i signalSet[i]->dimension()
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

	jointSignal:
	The joint signal must contain the same data 
	as the marginal signals, arranged subsequently
	in subdimension ranges. For example, 
	dimension 0 contains 1-d signalSet[0],
	dimensions 1-2 contain 2-d signalSet[1], and
	dimension 3 contains 1-d signalSet[2].
	
	signalSet:
	A set of marginal signals. They must all
	have the same size. However, their dimension
	can vary.

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

	TIMCORE real mutualInformation(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest = 1,
		real maxRelativeError = 0);

	TIMCORE real mutualInformation2(
		const std::vector<SignalPtr>& signalSet,
		integer kNearest = 1,
		real maxRelativeError = 0);

	//! Computes mutual information between 1-d marginal signals.
	/*!
	This is a convenience function that calls the
	more general mutualInformation() by splitting
	1-dimensional marginal signals from the joint signal.
	See the documentation for the more general mutualInformation().
	*/

	TIMCORE real mutualInformation(
		const SignalPtr& jointSignal,
		integer kNearest = 1,
		real maxRelativeError = 0);

}

#endif
