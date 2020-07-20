// Description: Kullback-Leibler divergence estimation
// Detail: Wang-Kulkarni-Verdu nearest neighbor estimator

#ifndef TIM_DIVERGENCE_WKV_H
#define TIM_DIVERGENCE_WKV_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include "pastel/sys/range.h"

namespace Tim
{

	//! Computes Kullback-Leibler divergence between signals.
	/*!
	xSignalSet:
	A set of signals representing trials for X.

	ySignalSet:
	A set of signals representing trials for X.

	returns:
	The Kullback-Leibler divergence between the signals.
	If the estimate is undefined, a NaN is returned.
	*/
	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	dreal divergenceWkv(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet);

	//! Computes Kullback-Leibler divergence between signals.
	/*!
	This is a convenience function that calls the
	more general divergenceWkv().

	See the documentation for that function.
	*/
	TIM dreal divergenceWkv(
		const Signal& xSignal,
		const Signal& ySignal);

}

#include "tim/core/divergence_wkv.hpp"

#endif
