#ifndef TIM_SIGNAL_PROPERTIES_H
#define TIM_SIGNAL_PROPERTIES_H

#include <pastel/sys/forwardrange.h>

namespace Tim
{

	template <typename SignalPtr_Iterator, typename Integer_Iterator>
	Integer2 sharedTimeInterval(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet);

	//! Returns the minimum number of samples among the signals.
	/*!
	Preconditions:
	SignalPtr_Iterator must dereference to SignalPtr.
	*/
	template <typename SignalPtr_Iterator>
	integer minSamples(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same dimension.
	/*!
	Preconditions:
	SignalPtr_Iterator must dereference to SignalPtr.
	*/
	template <typename SignalPtr_Iterator>
	bool equalDimension(
		const ForwardRange<SignalPtr_Iterator>& signalSet);

}

#include "tim/core/signal_properties.hpp"

#endif
