#ifndef TIM_TRANSFER_ENTROPY_HPP
#define TIM_TRANSFER_ENTROPY_HPP

#include "tim/core/transfer_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/null_iterator.h>

namespace Tim
{

	namespace Detail_TransferEntropy
	{

		template <
			typename X_Signal_Range,
			typename Y_Signal_Range,
			typename W_Signal_Range,
			typename Real_Filter_Iterator>
		real transferEntropy(
			const X_Signal_Range& xSignalSet,
			const Y_Signal_Range& ySignalSet,
			const W_Signal_Range& wSignalSet,
			integer timeWindowRadius,
			Signal* result,
			integer xLag, integer yLag, integer wLag,
			integer kNearest,
			const boost::iterator_range<Real_Filter_Iterator>& filter)
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

			std::vector<Signal> jointSignalSet;
			jointSignalSet.reserve(trials);

			Array<Signal> signalSet(Vector2i(trials, 3));
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
		typename W_Signal_Range,
		typename Real_Filter_Iterator>
	Signal temporalTransferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const W_Signal_Range& wSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer wLag,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
	{
		Signal result;
		Tim::Detail_TransferEntropy::transferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			timeWindowRadius,
			&result,
			xLag, yLag, wLag,
			kNearest,
			filter);
		return result;
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename W_Signal_Range>
	Signal temporalTransferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const W_Signal_Range& wSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer wLag,
		integer kNearest)
	{
		return Tim::temporalTransferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			timeWindowRadius,
			xLag, yLag, wLag,
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename W_Signal_Range>
	real transferEntropy(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const W_Signal_Range& wSignalSet,
		integer xLag, integer yLag, integer wLag,
		integer kNearest)
	{
		return Tim::Detail_TransferEntropy::transferEntropy(
			xSignalSet, ySignalSet, wSignalSet,
			0, 0,
			xLag, yLag, wLag,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
