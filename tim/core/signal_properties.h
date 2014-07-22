// Description: Properties of signals

#ifndef TIM_SIGNAL_PROPERTIES_H
#define TIM_SIGNAL_PROPERTIES_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/sys/tuple.h>

namespace Tim
{

	//! Returns the largest time interval on which all signals are defined.
	template <typename SignalPtr_Range, typename Integer_Iterator>
	Integer2 sharedTimeInterval(
		const SignalPtr_Range& signalSet,
		const boost::iterator_range<Integer_Iterator>& lagSet);

	//! Returns the largest time interval on which all signals are defined.
	/*!
	This is a convenience function that calls
	
	sharedTimeInterval(
		signalSet,
		constantRange(0, signalSet.size()));
	*/
	template <typename SignalPtr_Range>
	Integer2 sharedTimeInterval(
		const SignalPtr_Range& signalSet);

	//! Returns the minimum number of samples among the signals.
	template <typename SignalPtr_Range>
	integer minSamples(
		const SignalPtr_Range& signalSet);

	//! Returns true if all signals have the same dimension.
	template <typename SignalPtr_Range>
	bool equalDimension(
		const SignalPtr_Range& signalSet);

	//! Returns true if all signals have the same number of samples.
	template <typename SignalPtr_Range>
	bool equalSamples(
		const SignalPtr_Range& signalSet);

}

#include "tim/core/signal_properties.hpp"

#endif
