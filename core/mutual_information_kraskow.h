// Description: Mutual information estimation via entropy combination

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
	SignalPtr_X_Iterator dereferences to SignalPtr.
	SignalPtr_Y_Iterator dereferences to SignalPtr.
	Real_OutputIterator dereferences to a convertible to real.
	yLag >= 0
	timeWindowRadius >= 0
	kNearest > 0

	xSignalSet:
	A set of measurements (trials) of signal A. 

	ySignalSet:
	A set of measurements (trials) of signal B. 

	yLag:
	The delay in samples that is applied to signal B.

	timeWindowRadius:
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
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename Real_OutputIterator>
	void temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag = 0,
		integer kNearest = 1);

	template <typename Real_OutputIterator>
	void temporalMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag = 0,
		integer kNearest = 1);

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real mutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer yLag = 0,
		integer kNearest = 1);

	TIMCORE real mutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer yLag = 0,
		integer kNearest = 1);

}

#include "tim/core/mutual_information_kraskow.hpp"

#endif
