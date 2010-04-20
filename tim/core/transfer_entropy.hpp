#ifndef TIM_TRANSFER_ENTROPY_HPP
#define TIM_TRANSFER_ENTROPY_HPP

#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constantiterator.h>
#include <pastel/sys/nulliterator.h>

namespace Tim
{

	namespace Detail_TransferEntropy
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename SignalPtr_W_Iterator,
			typename Real_OutputIterator,
			typename Real_Filter_Iterator>
		real transferEntropy(
			const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
			const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
			const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
			integer timeWindowRadius,
			Real_OutputIterator result,
			integer xLag, integer yLag, integer wLag,
			integer kNearest,
			const ForwardRange<Real_Filter_Iterator>& filter,
			bool wantTemporal)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			ENSURE(odd(filter.size()));
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE_OP(xSignalSet.size(), ==, wSignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));

			if (xSignalSet.empty())
			{
				return 0;
			}

			const integer trials = xSignalSet.size();

			// Form the joint signal. Note the signals 
			// are merged in wXY order.

			std::vector<SignalPtr> jointSignalSet;
			jointSignalSet.reserve(trials);

			Array<SignalPtr> signalSet(trials, 3);
			std::copy(wSignalSet.begin(), wSignalSet.end(), signalSet.rowBegin(0));
			std::copy(xSignalSet.begin(), xSignalSet.end(), signalSet.rowBegin(1));
			std::copy(ySignalSet.begin(), ySignalSet.end(), signalSet.rowBegin(2));

			integer lagSet[] = {wLag, xLag, yLag};

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
		typename SignalPtr_W_Iterator,
		typename Real_OutputIterator,
		typename Real_Filter_Iterator>
	integer temporalTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag, integer yLag, integer wLag,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter)
	{
		return Tim::Detail_TransferEntropy::transferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			timeWindowRadius,
			result,
			xLag, yLag, wLag,
			kNearest,
			filter,
			true);
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_W_Iterator,
		typename Real_OutputIterator>
	integer temporalTransferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag, integer yLag, integer wLag,
		integer kNearest)
	{
		return Tim::temporalTransferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			timeWindowRadius,
			result,
			xLag, yLag, wLag,
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_W_Iterator>
	real transferEntropy(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		const ForwardRange<SignalPtr_W_Iterator>& wSignalSet,
		integer xLag, integer yLag, integer wLag,
		integer kNearest)
	{
		return Tim::Detail_TransferEntropy::transferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			0,
			NullIterator(),
			xLag, yLag, wLag,
			kNearest,
			constantRange((real)1, 1),
			false);
	}

}

#endif
