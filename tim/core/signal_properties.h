// Description: Properties of signals

#ifndef TIM_SIGNAL_PROPERTIES_H
#define TIM_SIGNAL_PROPERTIES_H

#include "tim/core/mytypes.h"
#include "tim/core/signal.h"

#include <pastel/sys/range.h>
#include <pastel/sys/tuple.h>

#include <algorithm>

namespace Tim
{

	//! Returns the largest time interval on which all signals are defined.
	template <
		ranges::forward_range Signal_Range, 
		ranges::forward_range Lag_Range>
	Integer2 sharedTimeInterval(
		const Signal_Range& signalSet,
		const Lag_Range& lagSet)
	{
		// In the following we think of having signals
		// embedded on the time axis, delayed with the given
		// lags. We accept to the merged signal only that
		// part in which all signals are present. E.g,
		// if the following are the time intervals that
		// three signals span:
		// 
		//        +--------------+
		//   +--------+
		//         +------+
		//
		// Then the their merged signal spans the following
		// time interval:
		//
		//         +--+
		//
		// Note that each signal also has its own lag which
		// must be added to the given lag.

		// Find out the common time interval
		// [tLeftMax, tRightMin[.

		auto lagIter = std::begin(lagSet);
		auto lagIterEnd = std::end(lagSet);
		auto signalIter = ranges::begin(signalSet);

		integer tLeftMax = (*lagIter) + signalIter->t();
		integer tRightMin = tLeftMax + signalIter->samples();
		++lagIter;
		++signalIter;

		while(lagIter != lagIterEnd)
		{
			const integer tLeft = *lagIter + signalIter->t();
			const integer tRight = tLeft + signalIter->samples();

			if (tLeft > tLeftMax)
			{
				tLeftMax = tLeft;
			}
			if (tRight < tRightMin)
			{
				tRightMin = tRight;
			}

			++lagIter;
			++signalIter;
		}

		if (tRightMin < tLeftMax)
		{
			return Integer2(0, 0);
		}
		
		return Integer2(tLeftMax, tRightMin);
	}

	//! Returns the largest time interval on which all signals are defined.
	/*!
	This is a convenience function that calls
	
	sharedTimeInterval(
		signalSet,
		constantRange(0, ranges::size(signalSet)));
	*/
	template <ranges::forward_range Signal_Range>
	Integer2 sharedTimeInterval(const Signal_Range& signalSet)
	{
		return sharedTimeInterval(
			signalSet,
			constantRange(0, ranges::size(signalSet)));
	}

	//! Returns the minimum number of samples among the signals.
	template <ranges::forward_range Signal_Range>
	integer minSamples(const Signal_Range& signalSet)
	{
		return ranges::empty(signalSet) ? 0 : ranges::min(
			ranges::views::transform(signalSet, 
			[](auto&& signal) {return signal.samples();}));
	}

	//! Returns true if all signals have the same dimension.
	template <ranges::forward_range Signal_Range>
	bool equalDimension(const Signal_Range& signalSet)
	{
		if (ranges::empty(signalSet))
		{
			return true;
		}

		integer first = std::begin(signalSet)->dimension();
		return ranges::all_of(
			ranges::views::transform(signalSet, 
			[](auto&& signal) {return signal.dimension();}), 
			[first](int x) {return x == first;});
	}

	//! Returns true if all signals have the same number of samples.
	template <ranges::forward_range Signal_Range>
	bool equalSamples(const Signal_Range& signalSet)
	{
		if (ranges::empty(signalSet))
		{
			return true;
		}

		integer first = std::begin(signalSet)->samples();
		return ranges::all_of(
			ranges::views::transform(signalSet, 
			[](auto&& signal) {return signal.samples();}), 
			[first](int x) {return x == first;});
	}

}

#endif
