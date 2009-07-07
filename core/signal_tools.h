#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"

#include <pastel/math/matrix.h>
#include <pastel/math/cholesky_decomposition.h>

#include <pastel/sys/smallset.h>

#include <iostream>

namespace Tim
{

	TIMCORE std::ostream& operator<<(std::ostream& stream, const Signal& signal);

	//! Generates uniform random variables in [0, 1]^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/

	TIMCORE SignalPtr generateUniform(
		integer samples,
		integer dimension);

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0
	*/

	TIMCORE SignalPtr generateGaussian(
		integer samples,
		integer dimension);

	//! Generates correlated gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0

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
		integer samples,
		integer dimension,
		const CholeskyDecompositionD& covarianceCholesky);

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale);

	TIMCORE SignalPtr mergeMarginal(
		const std::vector<SignalPtr>& signalList);

	TIMCORE void splitMarginal(
		const SignalPtr& jointSignal,
		std::vector<SignalPtr>& marginalSet);

	TIMCORE void splitMarginal(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		std::vector<SignalPtr>& marginalSet);

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& resultPointSet);

	TIMCORE void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result);

	TIMCORE void normalizeCovariance(
		const SignalPtr& signal,
		const MatrixD& covariance);

}

#endif
