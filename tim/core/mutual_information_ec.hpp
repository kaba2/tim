#ifndef TIM_MUTUAL_INFORMATION_EC_HPP
#define TIM_MUTUAL_INFORMATION_EC_HPP

#include "tim/core/mutual_information_ec.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/null_iterator.h>

namespace Tim
{

	namespace Detail_MutualInformation
	{

		template <
			typename Signal_X_Iterator,
			typename Signal_Y_Iterator,
			typename Real_Filter_Iterator>
		real mutualInformation(
			const boost::iterator_range<Signal_X_Iterator>& xSignalSet,
			const boost::iterator_range<Signal_Y_Iterator>& ySignalSet,
			integer timeWindowRadius,
			Signal* result,
			integer xLag, integer yLag,
			integer kNearest,
			const boost::iterator_range<Real_Filter_Iterator>& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			PENSURE_OP(xSignalSet.size(), ==, ySignalSet.size());
			PENSURE(equalDimension(xSignalSet));
			PENSURE(equalDimension(ySignalSet));
			ENSURE(odd(filter.size()));

			if (xSignalSet.empty() || ySignalSet.empty())
			{
				return 0;
			}

			// Copy the signals in an array.

			const integer trials = xSignalSet.size();

			Array<Signal> signalSet(Vector2i(trials, 2));
			std::copy(xSignalSet.begin(), xSignalSet.end(),
				signalSet.rowBegin(0));
			std::copy(ySignalSet.begin(), ySignalSet.end(),
				signalSet.rowBegin(1));

			// Describe the marginal signals.

			const Integer3 rangeSet[] = 
			{
				Integer3(0, 1, 1),
				Integer3(1, 2, 1)
			};

			const integer lagSet[] = {xLag, yLag};

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
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator,
		typename Real_Filter_Iterator>
	Signal temporalMutualInformation(
		const boost::iterator_range<Signal_X_Iterator>& xSignalSet,
		const boost::iterator_range<Signal_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
	{
		Signal result;
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
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator>
	Signal temporalMutualInformation(
		const boost::iterator_range<Signal_X_Iterator>& xSignalSet,
		const boost::iterator_range<Signal_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::temporalMutualInformation(
			xSignalSet, ySignalSet,
			timeWindowRadius,
			xLag, yLag,
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename Signal_X_Iterator,
		typename Signal_Y_Iterator>
	real mutualInformation(
		const boost::iterator_range<Signal_X_Iterator>& xSignalSet,
		const boost::iterator_range<Signal_Y_Iterator>& ySignalSet,
		integer xLag, integer yLag,
		integer kNearest)
	{
		return Tim::Detail_MutualInformation::mutualInformation(
			xSignalSet, ySignalSet,
			0,
			0,
			xLag, yLag,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
