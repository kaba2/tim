#ifndef TIM_MUTUAL_INFORMATION_KRASKOW_HPP
#define TIM_MUTUAL_INFORMATION_KRASKOW_HPP

#include "tim/core/mutual_information_kraskow.h"
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
			typename Signal_X_Iterator,
			typename Signal_Y_Iterator,
			typename Real_OutputIterator>
		real mutualInformation(
			const ForwardRange<Signal_X_Iterator>& xSignalSet,
			const ForwardRange<Signal_Y_Iterator>& ySignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer yLag,
			integer kNearest,
			bool wantTemporal)
		{
			ENSURE_OP(yLag, >=, 0);
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			const integer trials = xSignalSet.size();

			// Form the joint signal.

			std::vector<SignalPtr> jointSignalSet;
			jointSignalSet.reserve(trials);

			merge(xSignalSet, ySignalSet, 
				std::back_inserter(jointSignalSet), yLag);

			// Describe the marginal signals.

			const integer xBegin = 0;
			const integer xEnd = xSignalSet.front()->dimension();
			const integer yBegin = xEnd;
			const integer yEnd = yBegin + ySignalSet.front()->dimension();

			Integer3 rangeSet[2] = 
			{
				Integer3(xBegin, xEnd, 1),
				Integer3(yBegin, yEnd, 1)
			};

			if (wantTemporal)
			{
				temporalEntropyCombination(
					forwardRange(jointSignalSet.begin(), jointSignalSet.end()),
					forwardRange(rangeSet),
					timeWindowRadius,
					result,
					kNearest);
				
				return 0;
			}

			return entropyCombination(
				forwardRange(jointSignalSet.begin(), jointSignalSet.end()),
				forwardRange(rangeSet),
				kNearest);
		}

	}

	template <
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator,
		typename Real_OutputIterator>
	void temporalMutualInformation(
		const ForwardRange<Signal_X_Iterator>& xSignalSet,
		const ForwardRange<Signal_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag,
		integer kNearest)
	{
		Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet,
			ySignalSet,
			timeWindowRadius,
			result,
			yLag,
			kNearest,
			true);
	}

	template <typename Real_OutputIterator>
	void temporalMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag,
		integer kNearest)
	{
		Tim::temporalMutualInformation(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			timeWindowRadius,
			result,
			yLag,
			kNearest);
	}

	template <
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator>
	real mutualInformation(
		const ForwardRange<Signal_X_Iterator>& xSignalSet,
		const ForwardRange<Signal_Y_Iterator>& ySignalSet,
		integer yLag,
		integer kNearest)
	{
		return Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet,
			ySignalSet,
			0,
			NullIterator(),
			yLag,
			kNearest,
			false);
	}

}

#endif
