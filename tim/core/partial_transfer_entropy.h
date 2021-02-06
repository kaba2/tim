// Description: Partial transfer entropy estimation.

#ifndef TIM_PARTIAL_TRANSFER_ENTROPY_H
#define TIM_PARTIAL_TRANSFER_ENTROPY_H

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/range.h>
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

	//! Computes temporal partial transfer entropy.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)
	ranges::size(zSignalSet) == ranges::size(xSignalSet)
	wSignalSet.size() == wSignalSet.size()

	xSignalSet, ySignalSet, zSignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, Z, and W, respectively.

	xLag, yLag, zLag, wLag:
	The delays in samples that are applied to
	signals X, Y, Z, and W, respectively.

	timeWindowRadius:
	The radius of a time-window in samples over which
	the partial mutual information is estimated for a given
	time instant. Smaller values give sensitivity to
	temporal changes in partial mutual information, while larger 
	values give smaller variance.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range,
		typename W_Signal_Range,
		ranges::forward_range Filter_Range>
	Signal temporalPartialTransferEntropy(
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

	//! Computes temporal partial mutual information.
	/*!
	This is a convenience function that calls:

	temporalPartialTransferEntropy(
		xSignalSet, ySignalSet, zSignalSet, wSignalSet,
		timeWindowRadius,
		xLag, yLag, zLag, wLag,
		kNearest,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
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
		integer xLag = 0, integer yLag = 0,
		integer zLag = 0, integer wLag = 0,
		integer kNearest = 1)
	{
		return Tim::temporalPartialTransferEntropy(
			xSignalSet, ySignalSet,
			zSignalSet, wSignalSet,
			timeWindowRadius,
			xLag, yLag, zLag, wLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

	//! Computes partial transfer entropy.
	/*!
	Preconditions:
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)
	ranges::size(zSignalSet) == ranges::size(xSignalSet)
	wSignalSet.size() == ranges::size(xSignalSet)

	xSignalSet, ySignalSet, zSignalSet, wSignalSet:
	A set of measurements (trials) of signals
	X, Y, Z, and W, respectively.

	xLag, yLag, zLag, wLag:
	The delays in samples that are applied to
	signals X, Y, Z, and W, respectively.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

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
		integer xLag = 0, integer yLag = 0,
		integer zLag = 0, integer wLag = 0,
		integer kNearest = 1)
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
