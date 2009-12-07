// Description: Divergence estimation
// Detail: Wang-Kulkarni-Verdu nearest neighbor estimator

#ifndef PASTEL_DIVERGENCE_WKV_H
#define PASTEL_DIVERGENCE_WKV_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include "pastel/sys/forwardrange.h"

namespace Tim
{

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real divergenceWkv(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet);

	TIM real divergenceWkv(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal);

}

#include "tim/core/divergence_wkv.hpp"

#endif
