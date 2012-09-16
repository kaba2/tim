#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_HPP
#define TIM_PARTIAL_MUTUAL_INFORMATION_HPP

#include "tim/core/partial_mutual_information.h"
#include "tim/core/signal_tools.h"
#include "tim/core/signalpointset.h"
#include "tim/core/entropy_combination.h"
#include "tim/core/entropy_combination_t.h"

#include <pastel/sys/constant_iterator.h>
#include <pastel/sys/null_iterator.h>

namespace Tim
{

	namespace Detail_PartialMutualInformation
	{

		template <
			typename SignalPtr_X_Iterator,
			typename SignalPtr_Y_Iterator,
			typename SignalPtr_Z_Iterator,
			typename Real_Filter_Iterator>
		real partialMutualInformation(
			const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
			const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet,
			const boost::iterator_range<SignalPtr_Z_Iterator>& zSignalSet,
			integer timeWindowRadius,
			SignalPtr* result,
			integer xLag, integer yLag, integer zLag,
			integer kNearest,
			const boost::iterator_range<Real_Filter_Iterator>& filter)
		{
			ENSURE_OP(timeWindowRadius, >=, 0);
			ENSURE_OP(kNearest, >, 0);
			ENSURE(odd(filter.size()));
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

			// Note the signals are listed in XZY order.

			Array<SignalPtr> signalSet(Vector2i(trials, 3));
			std::copy(xSignalSet.begin(), xSignalSet.end(), signalSet.rowBegin(0));
			std::copy(zSignalSet.begin(), zSignalSet.end(), signalSet.rowBegin(1));
			std::copy(ySignalSet.begin(), ySignalSet.end(), signalSet.rowBegin(2));

			const integer lagSet[] = {xLag, zLag, yLag};

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
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator,
		typename Real_Filter_Iterator>
	SignalPtr temporalPartialMutualInformation(
		const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
		const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet,
		const boost::iterator_range<SignalPtr_Z_Iterator>& zSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag,
		integer kNearest,
		const boost::iterator_range<Real_Filter_Iterator>& filter)
	{
		SignalPtr result;
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
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator>
	SignalPtr temporalPartialMutualInformation(
		const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
		const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet,
		const boost::iterator_range<SignalPtr_Z_Iterator>& zSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag,
		integer kNearest)
	{
		return temporalPartialMutualInformation(
			xSignalSet, ySignalSet, zSignalSet, 
			timeWindowRadius,
			xLag, yLag, zLag, 
			kNearest,
			constantRange((real)1, 1));
	}

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename SignalPtr_Z_Iterator>
	real partialMutualInformation(
		const boost::iterator_range<SignalPtr_X_Iterator>& xSignalSet,
		const boost::iterator_range<SignalPtr_Y_Iterator>& ySignalSet,
		const boost::iterator_range<SignalPtr_Z_Iterator>& zSignalSet,
		integer xLag, integer yLag, integer zLag,
		integer kNearest)
	{
		return Tim::Detail_PartialMutualInformation::partialMutualInformation(
			xSignalSet, ySignalSet, zSignalSet,
			0, 0,
			xLag, yLag, zLag,
			kNearest,
			constantRange((real)1, 1));
	}

}

#endif
