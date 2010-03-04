// Description: Partial transfer entropy estimation.

#ifndef TIM_PARTIAL_TRANSFER_ENTROPY_H
#define TIM_PARTIAL_TRANSFER_ENTROPY_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes temporal partial transfer entropy.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()
	zSignalSet.size() == xSignalSet.size()
	wSignalSet.size() == wSignalSet.size()

	xSignalSet, ySignalSet, zSignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, Z, and W, respectively.

	xLag, yLag, zLag, wLag:
	The delays in samples that are applied to
	signals X, Y, Z, and W, respectively.

	timeWindowRadius:
	The radius of a time-window in samples over which
	the partial mutual information is estimated for a given
	time instant. Smaller values give sensitivity to
	temporal changes in partial mutual information, while larger 
	values give smaller variance.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	Returns:
	The number of time instants that had an
	undefined estimate. If not all estimates
	were undefined, they were reconstructed from 
	the defined estimates using interpolation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_W_Iterator,
		typename Real_OutputIterator>
	integer temporalPartialTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag = 0,
		integer yLag = 0,
		integer zLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes temporal partial transfer entropy.
	/*!
	This is a convenience function that calls:

	temporalPartialTransferEntropy(
		constantRange(xSignal), 
		constantRange(ySignal), 
		constantRange(zSignal), 
		constantRange(wSignal), 
		timeWindowRadius, result,
		xLag, yLag, zLag, wLag, kNearest);

	See the documentation for that function.
	*/

	template <typename Real_OutputIterator>
	integer temporalPartialTransferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		const SignalPtr& wSignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag = 0,
		integer yLag = 0,
		integer zLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes partial transfer entropy.
	/*!
	Preconditions:
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()
	zSignalSet.size() == xSignalSet.size()
	wSignalSet.size() == xSignalSet.size()

	xSignalSet, ySignalSet, zSignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, Z, and W, respectively.

	xLag, yLag, zLag, wLag:
	The delays in samples that are applied to
	signals X, Y, Z, and W, respectively.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_W_Iterator>
	real partialTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer xLag = 0,
		integer yLag = 0,
		integer zLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes partial transfer entropy.
	/*!
	This is a convenience function that calls:

	partialTransferEntropy(
		constantRange(xSignal), 
		constantRange(ySignal), 
		constantRange(zSignal), 
		constantRange(wSignal), 
		xLag, yLag, zLag, wLag, kNearest);

	See the documentation for that function.
	*/

	TIM real partialTransferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		const SignalPtr& wSignal,
		integer xLag = 0,
		integer yLag = 0,
		integer zLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

}

#include "tim/core/partial_transfer_entropy.hpp"

#endif
