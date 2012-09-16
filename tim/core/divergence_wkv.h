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
	An ensemble of signals representing trials
	of the same experiment X.

	ySignalSet:
	An ensemble of signals representing trials
	of the same experiment Y.

	Returns:
	The Kullback-Leibler divergence between the signals.
	If the estimate is undefined, a NaN is returned.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real divergenceWkv(
		const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
		const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet);

	//! Computes Kullback-Leibler divergence between signals.
	/*!
	This is a convenience function that calls the
	more general divergenceWkv().

	See the documentation for that function.
	*/
	TIM real divergenceWkv(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal);

}

#include "tim/core/divergence_wkv.hpp"

#endif
