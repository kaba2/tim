#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>

namespace Tim
{

	//! Generates uniform random variables in [0, 1]^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateUniform(
		integer dimension,
		integer size);

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateGaussian(
		integer dimension,
		integer size);

	//! Generates correlated gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0

	The correlated gaussian random variable is given by
	multiplying a standard gaussian random variable
	with the lower triangular part of the cholesky decomposition 
	of the correlation matrix.

	If the given correlation matrix turns out not to
	be numerically positive definite then
	the function call is equivalent to calling
	generateGaussian() (resulting in the 
	correlation matrix being identity).
	*/

	TIMCORE SignalPtr generateCorrelatedGaussian(
		integer dimension,
		integer size,
		const DynamicMatrix& correlation);

	/*!
	Preconditions:
	!signal.empty()
	*/

	TIMCORE void splitDimensions(
		const SignalPtr& signal,
		std::vector<SignalPtr>& signalSet);

	TIMCORE SignalPtr mergeDimensions(
		const std::vector<SignalPtr>& signalList);

}

#endif
