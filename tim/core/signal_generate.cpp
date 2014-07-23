#include "tim/core/signal_generate.h"

#include <pastel/math/matrix_tools.h>

#include <pastel/sys/random.h>

#include <boost/bind.hpp>

namespace Tim
{

	TIM Signal generateUniform(
		integer samples,
		integer dimension)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		Signal signal = Signal(samples, dimension);

		Matrix<real>::Iterator iter = signal.data().begin();
		Matrix<real>::Iterator iterEnd = signal.data().end();

		while(iter != iterEnd)
		{

			*iter = 2 * random<real>() - 1;
			++iter;
		}

		return signal;
	}

	TIM Signal generateGaussian(
		integer samples,
		integer dimension)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		Signal signal = Signal(samples, dimension);

		std::generate(signal.data().begin(), signal.data().end(),
			boost::bind(randomGaussian<real>));

		return signal;
	}

	TIM Signal generateCorrelatedGaussian(
		integer samples,
		integer dimension,
		const CholeskyDecomposition<real>& covarianceCholesky)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);
		ENSURE_OP(covarianceCholesky.lower().width(), ==, dimension);
		ENSURE(covarianceCholesky.succeeded());

		Signal correlatedGaussian = generateGaussian(samples, dimension);

		correlatedGaussian.data() *= transpose(covarianceCholesky.lower());

		return correlatedGaussian;
	}

	TIM Signal generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale)
	{
		ENSURE_OP(dimension, >, 0);
		ENSURE_OP(samples, >=, 0);

		Signal signal = Signal(samples, dimension);

		std::generate(signal.data().begin(), signal.data().end(),
			boost::bind(randomGeneralizedGaussian<real>, shape, scale));

		return signal;
	}

	TIM void generateTimeVaryingCoupling(
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

		xSignal.data().setSize(samples, 1);
		ySignal.data().setSize(samples, 1);
		zSignal.data().setSize(samples, 1);

		if (samples == 0)
		{
			return;
		}

		integer couplingStart = samples / 3;

		const integer couplingEnd = (samples * 2) / 3;
		integer couplingSamples = couplingEnd - couplingStart;
		real cyclesPerSample = 

			(2 * constantPi<real>()) / couplingSamples;

		Matrix<real>::Iterator xIter = xSignal.data().begin();
		Matrix<real>::Iterator yIter = ySignal.data().begin();
		Matrix<real>::Iterator zIter = zSignal.data().begin();

		for (integer i = 0;i < samples;++i)
		{
			real couplingYx = 0;
			real couplingZy = 0;
			if (i >= couplingStart && i < couplingEnd)
			{
				const real t = cyclesPerSample * (i - couplingStart);
				couplingYx = std::sin(t);
				couplingZy = std::cos(t);
			}

			real xPrevious = 0;
			real yPrevious = 0;
			real zPrevious = 0;
			if (i >= 1)
			{
				xPrevious = *(xIter - 1);
				yPrevious = *(yIter - 1);
				zPrevious = *(zIter - 1);
			}
			
			real xHistory = 0;
			if (i >= yxShift)
			{
				xHistory = *(xIter - yxShift);
			}

			real yHistory = 0;
			if (i >= zyShift)
			{
				yHistory = *(yIter - zyShift);
			}

			*xIter = 0.4 * xPrevious + 
				randomGaussian<real>();
			*yIter = 0.5 * yPrevious + 
				couplingYx * std::sin(xHistory) + 
				randomGaussian<real>();
			*zIter = 0.5 * zPrevious + 
				couplingZy * std::sin(yHistory) + 
				randomGaussian<real>();

			++xIter;
			++yIter;
			++zIter;
		}
	}

}
