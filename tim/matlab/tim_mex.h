// Description: Functions to ease interfacing with Matlab
// Documentation: tim_matlab_impl.txt

#ifndef TIM_TIM_MEX_H
#define TIM_TIM_MEX_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/matlab/matlab_argument.h>

namespace Tim
{

	//! Retrieves a reference to a signal.
	SignalPtr asSignal(const mxArray* signal);

	//! Retrieves the signals in a cell-array.
	/*!
	Returns:
	The number of signals in the array.

	Note: The signals are reported in Matlab's 
	linearized order.
	*/
	template <typename SignalPtr_Iterator>
	integer getSignals(
		const mxArray* input,
		SignalPtr_Iterator output);

	//! Retrieves an array of signals.
	void getSignalArray(
		const mxArray* signalSetArray, 
		Array<SignalPtr>& signalSet);
	
}

#include "tim/matlab/tim_mex.hpp"

#endif
