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
	template <typename SignalPtr_Iterator, typename Integer_Iterator>
	Integer2 sharedTimeInterval(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet,
		const boost::iterator_range<Integer_Iterator>& lagSet);

	//! Returns the largest time interval on which all signals are defined.
	/*!
	This is a convenience function that calls
	
	sharedTimeInterval(
		signalSet,
		constantRange(0, signalSet.size()));
	*/
	template <typename SignalPtr_Iterator>
	Integer2 sharedTimeInterval(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet);

	//! Returns the minimum number of samples among the signals.
	template <typename SignalPtr_Iterator>
	integer minSamples(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same dimension.
	template <typename SignalPtr_Iterator>
	bool equalDimension(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet);

	//! Returns true if all signals have the same number of samples.
	template <typename SignalPtr_Iterator>
	bool equalSamples(
		const boost::iterator_range<SignalPtr_Iterator>& signalSet);

}

#include "tim/core/signal_properties.hpp"

#endif
