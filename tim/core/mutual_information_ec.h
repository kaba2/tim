// Description: Mutual information estimation via entropy combination

#ifndef TIM_MUTUAL_INFORMATION_EC_H
#define TIM_MUTUAL_INFORMATION_EC_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/smallset.h>
#include <pastel/sys/forwardrange.h>

namespace Tim
{

	//! Computes temporal mutual information.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()

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

	Returns:
	The number of time instants that had an
	undefined estimate. If not all estimates
	were undefined, they were reconstructed from 
	the defined estimates using interpolation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename Real_OutputIterator,
		typename Real_Filter_Iterator>
	integer temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag, integer yLag,
		integer kNearest,
		const ForwardRange<Real_Filter_Iterator>& filter);

	//! Computes temporal mutual information.
	/*!
	This is a convenience function that calls:

	temporalMutualInformation(
		xSignalSet,
		ySignalSet,
		timeWindowRadius,
		result,
		xLag, yLag,
		filter,
		kNearest,
		constantRange((real)1, 1));

	See the documentation for that function.
	*/

	template <
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator,
		typename Real_OutputIterator>
	integer temporalMutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer timeWindowRadius,
		Real_OutputIterator result,
		integer xLag = 0, integer yLag = 0,
		integer kNearest = 1);

	//! Computes mutual information.
	/*!
	Preconditions:
	kNearest > 0
	ySignalSet.size() == xSignalSet.size()

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
		typename SignalPtr_X_Iterator,
		typename SignalPtr_Y_Iterator>
	real mutualInformation(
		const ForwardRange<SignalPtr_X_Iterator>& xSignalSet,
		const ForwardRange<SignalPtr_Y_Iterator>& ySignalSet,
		integer xLag = 0, integer yLag = 0,
		integer kNearest = 1);

}

#include "tim/core/mutual_information_ec.hpp"

#endif
