#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_HPP
#define TIM_PARTIAL_MUTUAL_INFORMATION_HPP

#include "tim/core/partial_mutual_information.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/nulliterator.h>

namespace Tim
{

	namespace Detail_PartialMutualInformation
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename SignalPtr_Z_Iterator,
			typename Real_OutputIterator>
		real partialMutualInformation(
			const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
			const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
			const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer xLag,
			integer yLag,
			integer zLag,
			integer kNearest,
			bool wantTemporal)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE_OP(xSignalSet.size(), ==, zSignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			PENSURE(equalDimension(zSignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			const integer trials = xSignalSet.size();

			// Note the signals are listed in XZY order.

			Array<SignalPtr, 2> signalSet(trials, 3);
			std::copy(xSignalSet.begin(), xSignalSet.end(), signalSet.rowBegin(0));
			std::copy(zSignalSet.begin(), zSignalSet.end(), signalSet.rowBegin(1));
			std::copy(ySignalSet.begin(), ySignalSet.end(), signalSet.rowBegin(2));

			const integer lagSet[] = {xLag, zLag, yLag};

			// Describe the marginal signals.

			Integer3 rangeSet[] = 
			{
				Integer3(0, 2, 1),
				Integer3(1, 3, 1),
				Integer3(1, 2, -1)
			};

			if (wantTemporal)
			{
				return temporalEntropyCombination(
					signalSet,
					forwardRange(rangeSet),
					timeWindowRadius,
					result,
					forwardRange(lagSet),
					kNearest);
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
		typename SignalPtr_Z_Iterator,
		typename Real_OutputIterator>
	integer temporalPartialMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		return Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet,
			ySignalSet,
			zSignalSet,
			timeWindowRadius,
			result,
			xLag,
			yLag,
			zLag,
			kNearest,
			true);
	}

	template <typename Real_OutputIterator>
	integer temporalPartialMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		return Tim::temporalPartialMutualInformation(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			forwardRange(constantIterator(zSignal)),
			timeWindowRadius,
			result,
			xLag,
			yLag,
			zLag,
			kNearest);
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator>
	real partialMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		integer xLag,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		return Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet,
			ySignalSet,
			zSignalSet,
			0,
			NullIterator(),
			xLag,
			yLag,
			zLag,
			kNearest,
			false);
	}

}

#endif
