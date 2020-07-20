// Description: Partial mutual information estimation.

#ifndef TIM_PARTIAL_MUTUAL_INFORMATION_H
#define TIM_PARTIAL_MUTUAL_INFORMATION_H

#include "tim/core/signal.h"

#include <pastel/sys/range.h>

namespace Tim
{

	//! Computes temporal partial mutual information.
	/*!
	Preconditions:
	timeWindowRadius >= 0
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)
	ranges::size(zSignalSet) == ranges::size(xSignalSet)

	xSignalSet, ySignalSet, zSignalSet:
	A set of measurements (trials) of signals
	X, Y, and Z, respectively.

	xLag, yLag, zLag:
	The delays in samples that are applied to
	signals X, Y, and Z, respectively.

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
		ranges::forward_range Filter_Range>
	Signal temporalPartialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer timeWindowRadius,
		integer xLag, integer yLag, integer zLag,
		integer kNearest,
		const Filter_Range& filter);

	//! Computes temporal partial mutual information.
	/*!
	This is a convenience function that calls:

	temporalPartialMutualInformation(
		xSignalSet, 
		ySignalSet, 
		zSignalSet, 
		timeWindowRadius,
		xLag, yLag, zLag, 
		kNearest,
		constantRange((dreal)1, 1));

	See the documentation for that function.
	*/

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range>
	Signal temporalPartialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer timeWindowRadius,
		integer xLag = 0, integer yLag = 0, integer zLag = 0,
		integer kNearest = 1);

	//! Computes partial mutual information.
	/*!
	Preconditions:
	kNearest > 0
	ranges::size(ySignalSet) == ranges::size(xSignalSet)
	ranges::size(zSignalSet) == ranges::size(xSignalSet)

	xSignalSet, ySignalSet, zSignalSet:
	A set of measurements (trials) of signals
	X, Y, and Z, respectively.

	xLag, yLag, zLag:
	The delays in samples that are applied to
	signals X, Y, and Z, respectively.

	kNearest:
	The number of nearest neighbors to use in the estimation.

	If the number of samples varies between trials, 
	then the minimum number of samples among the trials
	is used.
	*/

	template <
		typename X_Signal_Range,
		typename Y_Signal_Range,
		typename Z_Signal_Range>
	dreal partialMutualInformation(
		const X_Signal_Range& xSignalSet,
		const Y_Signal_Range& ySignalSet,
		const Z_Signal_Range& zSignalSet,
		integer xLag = 0, integer yLag = 0, integer zLag = 0,
		integer kNearest = 1);

}

#include "tim/core/partial_mutual_information.hpp"

#endif
