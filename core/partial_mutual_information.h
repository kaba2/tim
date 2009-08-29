// Description: Partial mutual information estimation.

#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_H
#define TIM_PARTIAL_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes partial mutual information I(A, B | C).
	/*!
	Preconditions:
	ySignalSet.size() == xSignalSet.size()
	zSignalSet.size() == xSignalSet.size()
	kNearest > 0
	yLag >= 0
	zLag >= 0
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename Real_OutputIterator>
	void temporalPartialMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag = 0,
		integer zLag = 0,
		integer kNearest = 1);

	template <typename Real_OutputIterator>
	void temporalPartialMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag = 0,
		integer zLag = 0,
		integer kNearest = 1);

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator>
	real partialMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		integer yLag = 0,
		integer zLag = 0,
		integer kNearest = 1);

}

#include "tim/core/partial_mutual_information.hpp"

#endif
