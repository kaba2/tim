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
			typename Signal_X_Iterator,
			typename Signal_Y_Iterator,
			typename Signal_Z_Iterator,
			typename Real_OutputIterator>
		real partialMutualInformation(
			const ForwardRange<Signal_X_Iterator>& xSignalSet,
			const ForwardRange<Signal_Y_Iterator>& ySignalSet,
			const ForwardRange<Signal_Z_Iterator>& zSignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer yLag,
			integer zLag,
			integer kNearest,
			bool wantTemporal)
		{
			ENSURE_OP(yLag, >=, 0);
			ENSURE_OP(zLag, >=, 0);
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

			// Form the joint signal. Note the signals 
			// are merged in XZY order.

			std::vector<SignalPtr> jointSignalSet;
			jointSignalSet.reserve(trials);

			merge(xSignalSet, zSignalSet, ySignalSet,
				std::back_inserter(jointSignalSet), yLag, zLag);

			// Describe the marginal signals.

			const integer xBegin = 0;
			const integer xEnd = xSignalSet.front()->dimension();
			const integer zBegin = xEnd;
			const integer zEnd = zBegin + zSignalSet.front()->dimension();
			const integer yBegin = zEnd;
			const integer yEnd = yBegin + ySignalSet.front()->dimension();

			Integer3 rangeSet[3] = 
			{
				Integer3(xBegin, zEnd, 1),
				Integer3(zBegin, yEnd, 1),
				Integer3(zBegin, zEnd, -1)
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
		typename Signal_Z_Iterator,
		typename Real_OutputIterator>
	void temporalPartialMutualInformation(
		const ForwardRange<Signal_X_Iterator>& xSignalSet,
		const ForwardRange<Signal_Y_Iterator>& ySignalSet,
		const ForwardRange<Signal_Y_Iterator>& zSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet,
			ySignalSet,
			zSignalSet,
			timeWindowRadius,
			result,
			yLag,
			zLag,
			kNearest,
			true);
	}

	template <typename Real_OutputIterator>
	void temporalPartialMutualInformation(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		const SignalPtr& zSignal,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer yLag,
		integer zLag,
		integer kNearest)
	{
		Tim::temporalPartialMutualInformation(
			forwardRange(constantIterator(xSignal)),
			forwardRange(constantIterator(ySignal)),
			forwardRange(constantIterator(zSignal)),
			timeWindowRadius,
			result,
			yLag,
			zLag,
			kNearest);
	}

	template <
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator,
		typename Signal_Z_Iterator>
	real partialMutualInformation(
		const ForwardRange<Signal_X_Iterator>& xSignalSet,
		const ForwardRange<Signal_Y_Iterator>& ySignalSet,
		const ForwardRange<Signal_Z_Iterator>& zSignalSet,
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
			yLag,
			zLag,
			kNearest,
			false);
	}

}

#endif
