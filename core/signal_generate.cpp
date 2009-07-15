#include "tim/core/signal_generate.h"

#include <pastel/math/matrix_tools.h>

#include <pastel/sys/random_vector.h>

namespace Tim
{

	TIMCORE SignalPtr generateUniform(
		integer samples,
		integer dimension)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		MatrixD::Iterator iter = signal->data().begin();
		const MatrixD::Iterator iterEnd = signal->data().end();

		while(iter != iterEnd)
		{
			*iter = 2 * random<real>() - 1;
			++iter;
		}

		return signal;
	}

	TIMCORE SignalPtr generateGaussian(
		integer samples,
		integer dimension)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		MatrixD::Iterator iter = signal->data().begin();
		const MatrixD::Iterator iterEnd = signal->data().end();

		while(iter != iterEnd)
		{
			*iter = randomGaussian<real>();
			++iter;
		}

		return signal;
	}

	TIMCORE SignalPtr generateCorrelatedGaussian(
		integer samples,
		integer dimension,
		const CholeskyDecompositionD& covarianceCholesky)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);
		ENSURE_OP(covarianceCholesky.lower().width(), ==, dimension);
		ENSURE(covarianceCholesky.succeeded());

		SignalPtr correlatedGaussian = generateGaussian(samples, dimension);

		correlatedGaussian->data() *= transpose(covarianceCholesky.lower());

		return correlatedGaussian;
	}

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		MatrixD::Iterator iter = signal->data().begin();
		const MatrixD::Iterator iterEnd = signal->data().end();

		while(iter != iterEnd)
		{
			*iter = randomGeneralizedGaussian<real>(shape, scale);
			++iter;
		}

		return signal;
	}

}
