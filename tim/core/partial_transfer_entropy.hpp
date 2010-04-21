#ifndef TIM_PARTIAL_TRANSFER_ENTROPY_HPP
#define TIM_PARTIAL_TRANSFER_ENTROPY_HPP

#include "tim/core/partial_transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/nulliterator.h>

namespace Tim
{

	namespace Detail_PartialTransferEntropy
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename SignalPtr_Z_Iterator,
			typename SignalPtr_W_Iterator,
			typename Real_Filter_Iterator>
		real partialTransferEntropy(
			const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
			const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
			const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
			const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
			integer timeWindowRadius,
			SignalPtr* result,
			integer xLag, integer yLag,	integer zLag, integer wLag,
			integer kNearest,
			const ForwardRange<Real_Filter_Iterator>& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE_OP(xSignalSet.size(), ==, zSignalSet.size());
			PENSURE_OP(xSignalSet.size(), ==, wSignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			PENSURE(equalDimension(zSignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			const integer trials = xSignalSet.size();

			// Form the joint signal. Note the signals 
			// are merged in wXZY order.

			std::vector<SignalPtr> jointSignalSet;
			jointSignalSet.reserve(trials);

			Array<SignalPtr> signalSet(trials, 4);
			std::copy(wSignalSet.begin(), wSignalSet.end(), signalSet.rowBegin(0));
			std::copy(xSignalSet.begin(), xSignalSet.end(), signalSet.rowBegin(1));
			std::copy(zSignalSet.begin(), zSignalSet.end(), signalSet.rowBegin(2));
			std::copy(ySignalSet.begin(), ySignalSet.end(), signalSet.rowBegin(3));

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
					forwardRange(rangeSet),
					timeWindowRadius,
					forwardRange(lagSet),
					kNearest,
					filter);
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
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_W_Iterator,
		typename Real_Filter_Iterator>
	SignalPtr temporalPartialTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag, integer wLag,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter)
	{
		SignalPtr result;
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
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_W_Iterator>
	SignalPtr temporalPartialTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
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
			constantRange((real)1, 1));
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename SignalPtr_W_Iterator>
	real partialTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_Z_Iterator>& zSignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer xLag, integer yLag, integer zLag, integer wLag,
		integer kNearest)
	{
		return Tim::Detail_PartialTransferEntropy::partialTransferEntropy(
			xSignalSet, ySignalSet, zSignalSet, wSignalSet,
			0, 0,
			xLag, yLag, zLag, wLag,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
