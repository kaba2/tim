// Description: Generation of common signals

#ifndef TIM_SIGNAL_GENERATE_H
#define TIM_SIGNAL_GENERATE_H

#include "tim/core/signal.h"

#include <pastel/math/matrix/cholesky_decomposition.h>
#include <pastel/math/matrix/matrix.h>

#include <pastel/sys/random.h>

namespace Tim
{

	//! Generates uniform random variables in [0, 1]^n.
	/*!
	Preconditions:
	dimension > 0
	size >= 0
	*/
	inline TIM void generateUniform(Signal signal)
	{
		ranges::generate(signal.data().range(), [](){return 2 * random<dreal>() - 1;});
	}

	inline TIM SignalData generateUniform(integer dimension, integer samples)
	{
		SignalData signal(dimension, samples);
		generateUniform(signal);
		return signal;
	}

	//! Generates standard gaussian random variables in R^n.
	/*!
	Preconditions:
	dimension > 0
	samples >= 0
	*/
	inline TIM void generateGaussian(Signal signal)
	{
		ranges::generate(signal.data().range(), [](){return randomGaussian<dreal>();});
	}

	inline TIM SignalData generateGaussian(integer dimension, integer samples)
	{
		SignalData signal(dimension, samples);
		generateGaussian(signal);
		return signal;
	}

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
	inline TIM void generateCorrelatedGaussian(
		Signal signal,
		const CholeskyDecompositionInplace<dreal>& covarianceCholesky)
	{
		ENSURE_OP(covarianceCholesky.lower().cols(), ==, signal.dimension());
		ENSURE(covarianceCholesky.succeeded());

		generateGaussian(signal);

		asMatrix(signal.data()) *= asMatrix(covarianceCholesky.lower().transpose());
	}

	inline TIM SignalData generateCorrelatedGaussian(
		integer dimension, integer samples, 
		const CholeskyDecompositionInplace<dreal>& covarianceCholesky) {
		SignalData signal(dimension, samples);
		generateCorrelatedGaussian(signal, covarianceCholesky);
		return signal;
	}

	inline TIM void generateGeneralizedGaussian(
		Signal signal, dreal shape,	dreal scale)
	{
		ranges::generate(signal.data().range(), 
			[shape, scale](){return randomGeneralizedGaussian<dreal>(shape, scale);});
	}

	inline TIM SignalData generateGeneralizedGaussian(
		integer dimension, integer samples, dreal shape, dreal scale)
	{
		SignalData signal(dimension, samples);
		generateGeneralizedGaussian(signal, shape, scale);
		return signal;
	}

	//! Generates a signal with time-varying coupling.
	/*!
	Preconditions:
	samples >= 0
	yzShift >= 0
	zyShift >= 0

	The signals are divide into three time regions.
	In the first and the third time regions, there is
	no coupling between x, y, and z. However, in
	the second time region x drives y, and y drives z.
	The amplitudes of these drives are given by a sine
	wave for the x->y, and by a cosine wave for the
	y->z. Thus, those estimators which are sensitive
	to temporal changes in coupling (e.g. partial 
	transfer entropy) should give similar coupling curves.	
	*/
	inline TIM void generateTimeVaryingCoupling(
		integer samples,
		integer yxShift,
		integer zyShift,
		Signal& xSignal,
		Signal& ySignal,
		Signal& zSignal)
	{
		ENSURE_OP(samples, >=, 0);
		ENSURE_OP(yxShift, >=, 0);
		ENSURE_OP(zyShift, >=, 0);

		ENSURE_OP(xSignal.samples(), ==, ySignal.samples());
		ENSURE_OP(xSignal.samples(), ==, zSignal.samples());
		ENSURE_OP(xSignal.dimension(), ==, 1);
		ENSURE_OP(ySignal.dimension(), ==, 1);
		ENSURE_OP(zSignal.dimension(), ==, 1);

		if (xSignal.samples() == 0)
		{
			return;
		}

		integer couplingStart = samples / 3;

		const integer couplingEnd = (samples * 2) / 3;
		integer couplingSamples = couplingEnd - couplingStart;
		dreal cyclesPerSample = 

			(2 * constantPi<dreal>()) / couplingSamples;

		auto xi = std::begin(xSignal.data().range());
		auto yi = std::begin(ySignal.data().range());
		auto zi = std::begin(zSignal.data().range());

		for (integer i = 0;i < samples;++i)
		{
			dreal couplingYx = 0;
			dreal couplingZy = 0;
			if (i >= couplingStart && i < couplingEnd)
			{
				const dreal t = cyclesPerSample * (i - couplingStart);
				couplingYx = std::sin(t);
				couplingZy = std::cos(t);
			}

			dreal xPrevious = 0;
			dreal yPrevious = 0;
			dreal zPrevious = 0;
			if (i >= 1)
			{
				xPrevious = *(xi - 1);
				yPrevious = *(yi - 1);
				zPrevious = *(zi - 1);
			}
			
			dreal xHistory = 0;
			if (i >= yxShift)
			{
				xHistory = *(xi - yxShift);
			}

			dreal yHistory = 0;
			if (i >= zyShift)
			{
				yHistory = *(yi - zyShift);
			}

			*xi = 0.4 * xPrevious + 
				randomGaussian<dreal>();
			*yi = 0.5 * yPrevious + 
				couplingYx * std::sin(xHistory) + 
				randomGaussian<dreal>();
			*zi = 0.5 * zPrevious + 
				couplingZy * std::sin(yHistory) + 
				randomGaussian<dreal>();

			++xi;
			++yi;
			++zi;
		}
	}

}

#endif
