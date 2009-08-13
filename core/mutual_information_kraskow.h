// Description: Mutual information estimation via Kraskow's algorithm

#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_H
#define TIM_MUTUAL_INFORMATION_KRASKOW_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/smallset.h>
#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes mutual information between two signals.
	/*!
	Preconditions:
	Signal_A_Iterator dereferences to SignalPtr.
	Signal_B_Iterator dereferences to SignalPtr.
	Real_OutputIterator dereferences to a convertible to real.
	bLag >= 0
	sigma >= 0
	kNearest > 0
	
	aSignalSet:
	A set of measurements (trials) of signal A. 

	bSignalSet:
	A set of measurements (trials) of signal B. 

	bLag:
	The delay in samples that is applied to signal B.

	sigma:
	The radius of a time-window in samples over which
	the mutual information is estimated for a given
	time instant. Smaller values give sensitivity to
	temporal changes in mutual information, while larger 
	values give smaller variance.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used to compute the mi.

	If H(X) denotes the differential entropy of a real 
	random variable X, and we are given random variables
	X and Y, then mutual information is given by:
	
	I(X, Y) = H(X) + H(Y) - H(X, Y)

	See 'tim/core/differential_entropy.h' for more information
	on differential entropy.
	*/

	template <
		typename Signal_A_Iterator,
		typename Signal_B_Iterator,
		typename Real_OutputIterator>
	void mutualInformation(
		const ForwardRange<Signal_A_Iterator>& aSignalSet,
		const ForwardRange<Signal_B_Iterator>& bSignalSet,
		Real_OutputIterator result,
		integer bLag = 0,
		integer sigma = -1,
		integer kNearest = 1);

	template <
		typename Real_OutputIterator>
	void mutualInformation(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		Real_OutputIterator result,
		integer bLag = 0,
		integer sigma = -1,
		integer kNearest = 1);

}

#include "tim/core/mutual_information_kraskow.hpp"

#endif
