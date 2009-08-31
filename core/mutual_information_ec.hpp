#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_HPP
#define TIM_MUTUAL_INFORMATION_KRASKOW_HPP

#include "tim/core/mutual_information_ec.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/nulliterator.h>

namespace Tim
{

	namespace Detail_MutualInformation
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename Real_OutputIterator>
		real mutualInformation(
			const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
			const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer xLag,
			integer yLag,
			integer kNearest,
			bool wantTemporal)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			// Copy the signals in an array.

			const integer trials = xSignalSet.size();

			Array<SignalPtr, 2> signalSet(trials, 2);
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
				temporalEntropyCombination(
					signalSet, 
					forwardRange(rangeSet),
					timeWindowRadius,
					result,
					forwardRange(lagSet),
					kNearest);
				
				return 0;
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
		typename Real_OutputIterator>
	void temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag,
		integer yLag,
		integer kNearest)
	{
		Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet,
			ySignalSet,
			timeWindowRadius,
			result,
			xLag, yLag,
			kNearest,
			true);
	}

	template <typename Real_OutputIterator>
	void temporalMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag,
		integer yLag,
		integer kNearest)
	{
		Tim::temporalMutualInformation(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			timeWindowRadius,
			result,
			xLag, yLag,
			kNearest);
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
			xSignalSet,
			ySignalSet,
			0,
			NullIterator(),
			xLag, yLag,
			kNearest,
			false);
	}

}

#endif
