// Description: Functions to ease interfacing with Matlab
// Documentation: tim_matlab_impl.txt

#ifndef TIM_TIM_MEX_H
#define TIM_TIM_MEX_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/matlab/matlab_argument.h>

#include <vector>

namespace Tim
{

	//! Retrieves a reference to a signal.
	Signal asSignal(const mxArray* signal);

	//! Retrieves the signals in a cell-array.
	/*!
	Note: The signals are reported in Matlab's 
	linearized order.
	*/
	std::vector<Signal> getSignals(
		const mxArray* input);

	//! Retrieves an array of signals.
	Array<Signal> getSignalArray(
		const mxArray* signalSetArray);
	
}

#include "tim/corematlab/tim_mex.hpp"

#endif
