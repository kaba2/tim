// Description: Transfer entropy estimation.

#ifndef TIM_TRANSFER_ENTROPY_H
#define TIM_TRANSFER_ENTROPY_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes temporal transfer entropy.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()
	wSignalSet.size() == wSignalSet.size()

	xSignalSet, ySignalSet, zSignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, Z, and W, respectively.

	xLag, yLag, wLag:
	The delays in samples that are applied to
	signals X, Y, and W, respectively.

	timeWindowRadius:
	The radius of a time-window in samples over which
	the partial mutual information is estimated for a given
	time instant. Smaller values give sensitivity to
	temporal changes in partial mutual information, while larger 
	values give smaller variance.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_W_Iterator,
		typename Real_OutputIterator>
	void temporalTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag = 0,
		integer yLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes temporal transfer entropy.
	/*!
	This is a convenience function that calls:

	temporalTransferEntropy(
		forwardRange(constantIterator(xSignal)), 
		forwardRange(constantIterator(ySignal)), 
		forwardRange(constantIterator(wSignal)), 
		timeWindowRadius, result,
		xLag, yLag, wLag, kNearest);

	See the documentation for that function.
	*/

	template <typename Real_OutputIterator>
	void temporalTransferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& wSignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag = 0,
		integer yLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes transfer entropy.
	/*!
	Preconditions:
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()
	wSignalSet.size() == xSignalSet.size()

	xSignalSet, ySignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, and W, respectively.

	xLag, yLag, wLag:
	The delays in samples that are applied to
	signals X, Y, and W, respectively.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_W_Iterator>
	real transferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer xLag = 0,
		integer yLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

	//! Computes transfer entropy.
	/*!
	This is a convenience function that calls:

	transferEntropy(
		forwardRange(constantIterator(xSignal)), 
		forwardRange(constantIterator(ySignal)), 
		forwardRange(constantIterator(wSignal)), 
		xLag, yLag, wLag, kNearest);

	See the documentation for that function.
	*/

	TIM real transferEntropy(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& wSignal,
		integer xLag = 0,
		integer yLag = 0,
		integer wLag = 0,
		integer kNearest = 1);

}

#include "tim/core/transfer_entropy.hpp"

#endif
