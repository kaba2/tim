// Description: Mutual information estimation via entropy combination

#ifndef TIM_MUTUAL_INFORMATION_EC_H
#define TIM_MUTUAL_INFORMATION_EC_H

#include "tim/core/signal.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/math/matrix/matrix.h>
#include <pastel/math/matrix/cholesky_decomposition.h>
#include <pastel/sys/range.h>
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

	//! Computes temporal mutual information.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)

	xSignalSet, ySignalSet:
	A set of measurements (trials) of 
	signals X and Y, respectively. 

	xLag, yLag:
	The delays in samples that are applied to
	signals X and Y, respectively.

	timeWindowRadius:
	The radius of a time-window in samples over which
	the mutual information is estimated for a given
	time instant. Smaller values give sensitivity to
	temporal changes in mutual information, while larger 
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

	//! Computes temporal mutual information.
	/*!
	This is a convenience function that calls:

	temporalMutualInformation(
		xSignalSet,
		ySignalSet,
		timeWindowRadius,
		xLag, yLag,
		filter,
		kNearest,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/
	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	SignalData temporalMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		integer timeWindowRadius,
		integer xLag = 0, integer yLag = 0,
		integer kNearest = 1)
	{
		return Tim::temporalMutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			xLag, yLag,
			kNearest,
			constantRange((dreal)1, 1));
	}

	//! Computes mutual information.
	/*!
	Preconditions:
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)

	xSignalSet, ySignalSet:
	A set of measurements (trials) of 
	signals X and Y, respectively. 

	xLag, yLag:
	The delays in samples that are applied to
	signals X and Y, respectively.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/
	template <
		typename X_Signal_Range,
		typename Y_Signal_Range>
	dreal mutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		integer xLag = 0, integer yLag = 0,
		integer kNearest = 1)
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
