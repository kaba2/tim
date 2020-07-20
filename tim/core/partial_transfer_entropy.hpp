#ifndef TIM_PARTIAL_TRANSFER_ENTROPY_HPP
#define TIM_PARTIAL_TRANSFER_ENTROPY_HPP

#include "tim/core/partial_transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/iterator/null_iterator.h>

namespace Tim
{

	namespace Detail_PartialTransferEntropy
	{

		template <
			typename X_Signal_Range,
			typename Y_Signal_Range,
			typename Z_Signal_Range,
			typename W_Signal_Range,
			ranges::forward_range Filter_Range>
		dreal partialTransferEntropy(
			const X_Signal_Range& xSignalSet,
			const Y_Signal_Range& ySignalSet,
			const Z_Signal_Range& zSignalSet,
			const W_Signal_Range& wSignalSet,
			integer timeWindowRadius,
			Signal* result,
			integer xLag, integer yLag,	integer zLag, integer wLag,
			integer kNearest,
			const Filter_Range& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(ySignalSet));
			PENSURE_OP(ranges::size(xSignalSet), ==, ranges::size(zSignalSet));
			PENSURE_OP(ranges::size(xSignalSet), ==, wSignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			PENSURE(equalDimension(zSignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			integer trials = ranges::size(xSignalSet);

			// Form the joint signal. Note the signals 
			// are merged in wXZY order.

			std::vector<Signal> jointSignalSet;
			jointSignalSet.reserve(trials);

			Array<Signal> signalSet(Vector2i(trials, 4));
			std::copy(wSignalSet.begin(), wSignalSet.end(), signalSet.rowBegin(0));
			std::copy(std::begin(xSignalSet), std::end(xSignalSet), signalSet.rowBegin(1));
			std::copy(std::begin(zSignalSet), std::end(zSignalSet), signalSet.rowBegin(2));
			std::copy(std::begin(ySignalSet), std::end(ySignalSet), signalSet.rowBegin(3));

			integer lagSet[] = {wLag, xLag, zLag, yLag};

			// Describe the marginal signals.

			Integer3 rangeSet[] = 
			{
				Integer3(0, 3, 1),
				Integer3(1, 4, 1),
				Integer3(1, 3, -1)
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
		typename W_Signal_Range,
		ranges::forward_range Filter_Range>
	SignalData temporalPartialTransferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		const W_Signal_Range& wSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag, integer wLag,
		integer kNearest,
		const Filter_Range& filter)
	{
		SignalData result;
		Tim::Detail_PartialTransferEntropy::partialTransferEntropy(
			xSignalSet, ySignalSet, zSignalSet, wSignalSet,
			timeWindowRadius,
			&result,
			xLag, yLag, zLag, wLag,
			kNearest,
			filter);
		return result;
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range,
		typename W_Signal_Range>
	Signal temporalPartialTransferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		const W_Signal_Range& wSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag, integer wLag,
		integer kNearest)
	{
		return Tim::temporalPartialTransferEntropy(
			xSignalSet, ySignalSet,
			zSignalSet, wSignalSet,
			timeWindowRadius,
			xLag, yLag, zLag, wLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range,
		typename W_Signal_Range>
	dreal partialTransferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		const W_Signal_Range& wSignalSet,
		integer xLag, integer yLag, integer zLag, integer wLag,
		integer kNearest)
	{
		return Tim::Detail_PartialTransferEntropy::partialTransferEntropy(
			xSignalSet, ySignalSet, zSignalSet, wSignalSet,
			0, 0,
			xLag, yLag, zLag, wLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

}

#endif
