#ifndef TIM_MUTUAL_INFORMATION_EC_HPP
#define TIM_MUTUAL_INFORMATION_EC_HPP

#include "tim/core/mutual_information_ec.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/iterator/null_iterator.h>

namespace Tim
{

	namespace Detail_MutualInformation
	{

		template <
			typename X_Signal_Range,
			typename Y_Signal_Range,
			ranges::forward_range Filter_Range>
		dreal mutualInformation(
			const X_Signal_Range& xSignalSet,
			const Y_Signal_Range& ySignalSet,
			integer timeWindowRadius,
			SignalData* result,
			integer xLag, integer yLag,
			integer kNearest,
			const Filter_Range& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(ySignalSet));
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			ENSURE(odd(ranges::size(filter)));

			if (ranges::empty(xSignalSet) || ranges::empty(ySignalSet))
			{
				return 0;
			}

			// Copy the signals in an array.

			integer trials = ranges::size(xSignalSet);

			Array<Signal> signalSet(Vector2i(trials, 2));
			std::copy(std::begin(xSignalSet), std::end(xSignalSet),
				signalSet.rowBegin(0));
			std::copy(std::begin(ySignalSet), std::end(ySignalSet),
				signalSet.rowBegin(1));

			// Describe the marginal signals.

			Integer3 rangeSet[] = 
			{
				Integer3(0, 1, 1),
				Integer3(1, 2, 1)
			};

			integer lagSet[] = {xLag, yLag};

			// Compute entropy combination.

			if (result)
			{

				*result = temporalEntropyCombination(
					signalSet, 
					range(rangeSet),
					timeWindowRadius,
					range(lagSet),
					kNearest,
					filter);
				
				return 0;
			}

			return entropyCombination(
				signalSet,
				range(rangeSet),
				range(lagSet),
				kNearest);
		}

	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		ranges::forward_range Filter_Range>
	SignalData temporalMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag,
		integer kNearest,
		const Filter_Range& filter)
	{
		SignalData result;
		Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			&result,
			xLag, yLag,
			kNearest,
			filter);
		return result;
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	Signal temporalMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::temporalMutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			xLag, yLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	dreal mutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet, ySignalSet,
			0,
			0,
			xLag, yLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

}

#endif
