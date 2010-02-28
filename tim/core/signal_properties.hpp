#ifndef TIM_SIGNAL_PROPERTIES_HPP
#define TIM_SIGNAL_PROPERTIES_HPP

#include "tim/core/signal_properties.h"
#include "tim/core/mytypes.h"

#include <algorithm>

namespace Tim
{

	template <typename SignalPtr_Iterator, typename Integer_Iterator>
	Integer2 sharedTimeInterval(
		const ForwardRange<SignalPtr_Iterator>& signalSet,
		const ForwardRange<Integer_Iterator>& lagSet)
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

		// Find out the common time interval
		// [tLeftMax, tRightMin[.

		Integer_Iterator lagIter = lagSet.begin();
		const Integer_Iterator lagIterEnd = lagSet.end();
		SignalPtr_Iterator signalIter = signalSet.begin();

		integer tLeftMax = (*lagIter);
		integer tRightMin = tLeftMax + (*signalIter)->samples();
		++lagIter;
		++signalIter;

		while(lagIter != lagIterEnd)
		{
			const integer tLeft = *lagIter;
			const integer tRight = tLeft + (*signalIter)->samples();

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

	template <typename SignalPtr_Iterator>
	integer minSamples(
		const ForwardRange<SignalPtr_Iterator>& signalSet)
	{
		if (signalSet.empty())
		{
			return 0;
		}

		SignalPtr_Iterator iter = signalSet.begin();
		const SignalPtr_Iterator iterEnd = signalSet.end();

		integer samples = (*iter)->samples();
		++iter;

		while(iter != iterEnd)
		{
			samples = std::max(samples, (*iter)->samples());

			++iter;
		}

		return samples;
	}

	template <typename SignalPtr_Iterator>
	bool equalDimension(
		const ForwardRange<SignalPtr_Iterator>& signalSet)
	{
		if (signalSet.empty())
		{
			return true;
		}

		SignalPtr_Iterator iter = signalSet.begin();
		const SignalPtr_Iterator iterEnd = signalSet.end();

		integer dimension = signalSet.front()->dimension();
		++iter;

		while(iter != iterEnd)
		{
			if ((*iter)->dimension() != dimension)
			{
				return false;
			}

			++iter;
		}

		return true;
	}

}

#endif
