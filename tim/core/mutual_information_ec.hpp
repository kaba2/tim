#ifndef TIM_MUTUAL_INFORMATION_EC_HPP
#define TIM_MUTUAL_INFORMATION_EC_HPP

#include "tim/core/mutual_information_ec.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/nulliterator.h>

namespace Tim
{

	namespace Detail_MutualInformation
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename Real_OutputIterator,
			typename Real_Filter_Iterator>
		real mutualInformation(
			const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
			const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer xLag, integer yLag,
			integer kNearest,
			const ForwardRange<Real_Filter_Iterator>& filter,
			bool wantTemporal)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			ENSURE(odd(filter.size()));

			if (xSignalSet.empty() || ySignalSet.empty())
			{
				return 0;
			}

			// Copy the signals in an array.

			const integer trials = xSignalSet.size();

			Array<SignalPtr> signalSet(trials, 2);
			std::copy(xSignalSet.begin(), xSignalSet.end(),
				signalSet.rowBegin(0));
			std::copy(ySignalSet.begin(), ySignalSet.end(),
				signalSet.rowBegin(1));

			// Describe the marginal signals.

			const Integer3 rangeSet[] = 
			{
				Integer3(0, 1, 1),
				Integer3(1, 2, 1)
			};

			const integer lagSet[] = {xLag, yLag};

			// Compute entropy combination.

			if (wantTemporal)
			{
				return temporalEntropyCombination(
					signalSet, 
					forwardRange(rangeSet),
					timeWindowRadius,
					result,
					forwardRange(lagSet),
					kNearest,
					filter);
			}

			return entropyCombination(
				signalSet,
				forwardRange(rangeSet),
				forwardRange(lagSet),
				kNearest);
		}

	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename Real_OutputIterator,
		typename Real_Filter_Iterator>
	integer temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag,
		integer yLag,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter)
	{
		return Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			result,
			xLag, yLag,
			kNearest,
			filter,
			true);
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename Real_OutputIterator>
	integer temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::temporalMutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			result,
			xLag, yLag,
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real mutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet, ySignalSet,
			0,
			NullIterator(),
			xLag, yLag,
			kNearest,
			constantRange((real)1, 1),
			false);
	}

}

#endif
