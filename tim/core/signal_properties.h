#ifndef TIM_SIGNAL_PROPERTIES_H
#define TIM_SIGNAL_PROPERTIES_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/forwardrange.h>
#include <pastel/sys/tuple.h>

namespace Tim
{

	//! Counts the number of samples from the start which begin with a NaN.
	TIM integer nanSamples(const SignalPtr& signal);

	//! Returns the largest time interval on which all signals are defined.
	template <typename SignalPtr_Iterator, typename Integer_Iterator>
	Integer2 sharedTimeInterval(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet);

	//! Returns the largest time interval on which all signals are defined.
	/*!
	This is a convenience function that calls
	
	sharedTimeInterval(
		signalSet,
		constantRange(0, signalSet.size()));
	*/
	template <typename SignalPtr_Iterator>
	Integer2 sharedTimeInterval(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Returns the minimum number of samples among the signals.
	template <typename SignalPtr_Iterator>
	integer minSamples(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same dimension.
	template <typename SignalPtr_Iterator>
	bool equalDimension(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same number of samples.
	template <typename SignalPtr_Iterator>
	bool equalSamples(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

}

#include "tim/core/signal_properties.hpp"

#endif
