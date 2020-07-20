#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_HPP
#define TIM_PARTIAL_MUTUAL_INFORMATION_HPP

#include "tim/core/partial_mutual_information.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/iterator/null_iterator.h>

namespace Tim
{

	namespace Detail_PartialMutualInformation
	{

		template <
			typename X_Signal_Range,
			typename Y_Signal_Range,
			typename Z_Signal_Range,
			ranges::forward_range Filter_Range>
		dreal partialMutualInformation(
			const X_Signal_Range& xSignalSet,
			const Y_Signal_Range& ySignalSet,
			const Z_Signal_Range& zSignalSet,
			integer timeWindowRadius,
			Signal* result,
			integer xLag, integer yLag, integer zLag,
			integer kNearest,
			const Filter_Range& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			ENSURE(odd(ranges::size(filter)));
			PENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(ySignalSet));
			PENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(zSignalSet));
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			PENSURE(equalDimension(zSignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			integer trials = ranges::size(xSignalSet);

			// Note the signals are listed in XZY order.

			Array<Signal> signalSet(Vector2i(trials, 3));
			std::copy(std::begin(xSignalSet), std::end(xSignalSet), signalSet.rowBegin(0));
			std::copy(std::begin(zSignalSet), std::end(zSignalSet), signalSet.rowBegin(1));
			std::copy(std::begin(ySignalSet), std::end(ySignalSet), signalSet.rowBegin(2));

			integer lagSet[] = {xLag, zLag, yLag};

			// Describe the marginal signals.

			Integer3 rangeSet[] = 
			{
				Integer3(0, 2, 1),
				Integer3(1, 3, 1),
				Integer3(1, 2, -1)
			};

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
		typename Z_Signal_Range,
		ranges::forward_range Filter_Range>
	SignalData temporalPartialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag,
		integer kNearest,
		const Filter_Range& filter)
	{
		SignalData result;
		Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet, ySignalSet, zSignalSet,
			timeWindowRadius,
			&result,
			xLag, yLag,	zLag,
			kNearest,
			filter);
		return result;
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range>
	Signal temporalPartialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag,
		integer kNearest)
	{
		return temporalPartialMutualInformation(
			xSignalSet, ySignalSet, zSignalSet, 
			timeWindowRadius,
			xLag, yLag, zLag, 
			kNearest,
			constantRange((dreal)1, 1));
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range>
	dreal partialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer xLag, integer yLag, integer zLag,
		integer kNearest)
	{
		return Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet, ySignalSet, zSignalSet,
			0, 0,
			xLag, yLag, zLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

}

#endif
