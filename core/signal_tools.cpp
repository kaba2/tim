#include "tim/core/signal_tools.h"

#include <pastel/math/cholesky_decomposition.h>
#include <pastel/math/matrix_tools.h>

#include <pastel/sys/random_vector.h>
#include <pastel/sys/view_all.h>

namespace Tim
{

	TIMCORE std::ostream& operator<<(
		std::ostream& stream, const Signal& signal)
	{
		const integer dimension = signal.width();
		const integer samples = signal.height();
		if (dimension > 1)
		{
			for (integer i = 0;i < samples;++i)
			{
				for (integer j = 0;j < dimension;++j)
				{
					stream << signal[i][j] << " ";
				}
				stream << std::endl;
			}
		}
		else
		{
			for (integer i = 0;i < samples;++i)
			{
				stream << signal[i][0] << ", ";
			}
		}

		return stream;
	}

	TIMCORE SignalPtr generateUniform(
		integer samples,
		integer dimension)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		Signal::Iterator iter = signal->begin();
		const Signal::Iterator iterEnd = signal->end();

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
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		Signal::Iterator iter = signal->begin();
		const Signal::Iterator iterEnd = signal->end();

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
		const CholeskyDecompositionD& correlationCholesky)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);
		ENSURE(correlationCholesky.lower().width() == dimension);
		ENSURE(correlationCholesky.succeeded());

		SignalPtr correlatedGaussian = generateGaussian(samples, dimension);
	
		for (integer i = 0;i < samples;++i)
		{
			(*correlatedGaussian)[i] = 
				correlationCholesky.lower() * (*correlatedGaussian)[i];
		}

		return correlatedGaussian;
	}

	TIMCORE SignalPtr generateGeneralizedGaussian(
		integer samples,
		integer dimension,
		real shape,
		real scale)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);

		SignalPtr signal = SignalPtr(new Signal(samples, dimension));

		Signal::Iterator iter = signal->begin();
		const Signal::Iterator iterEnd = signal->end();

		while(iter != iterEnd)
		{
			*iter = randomGeneralizedGaussian<real>(shape, scale);
			++iter;
		}

		return signal;
	}

	TIMCORE void splitDimensions(
		const SignalPtr& jointSignal,
		std::vector<SignalPtr>& marginalSet)
	{
		const integer dimension = jointSignal->width();
		const integer n = jointSignal->height();

		for (integer x = 0;x < dimension;++x)
		{
			const SignalPtr signal = SignalPtr(new Signal(n, 1));
			for (integer y = 0;y < n;++y)
			{
				(*signal)(y) = (*jointSignal)(y, x);
			}

			marginalSet.push_back(signal);
		}
	}

	TIMCORE SignalPtr mergeDimensions(
		const std::vector<SignalPtr>& signalList)
	{
		if (signalList.empty())
		{
			return SignalPtr();
		}

		const integer size = signalList.front()->height();

		const integer signals = signalList.size();
		integer bigDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			bigDimension += signalList[i]->width();
			ENSURE2(signalList[0]->height() == signalList[i]->height(), 
				signalList[0]->height(), signalList[i]->height());
		}

		SignalPtr bigSignal(new Signal(size, bigDimension));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			copy(signalList[i]->constView(),
				subView(bigSignal->view(), 
				Rectangle2(0, dimensionOffset, size, 
				dimensionOffset + signalList[i]->width())));

			dimensionOffset += signalList[i]->width();
		}

		return bigSignal;
	}

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& resultPointSet)
	{
		const integer dimension = signal->width();
		const integer n = signal->height();
		
		std::vector<PointD> pointSet;
		pointSet.reserve(n);

		for (integer i = 0;i < n;++i)
		{
			pointSet.push_back(PointD((*signal)[i]));
		}

		pointSet.swap(resultPointSet);
	}

	TIMCORE void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result)
	{
		const integer dimension = signal->width();
		const integer samples = signal->height();

		result.setSize(dimension, dimension);
		result.set(0);

		const VectorD mean = sum(*signal) / samples;

		for (integer i = 0;i < samples;++i)
		{
			result += outerProduct((*signal)[i] - mean, (*signal)[i] - mean);
		}

		result /= samples;
	}

	TIMCORE void normalizeCovariance(
		const SignalPtr& signal,
		const MatrixD& covariance)
	{
		// Let X be the signal matrix with
		// each sample as a _column_.
		// Then the covariance C of the signal
		// is given by:
		//
		// C = X X^T
		//
		// Problem: find a matrix L
		// by which the signal X transforms
		// into a signal Y = LX 
		// having identity covariance.
		//
		// Solution:
		// 
		// Y Y^T = I
		// =>
		// (LX) (LX)^T = I
		// =>
		// L X X^T L^T = I
		// =>
		// L C L^T = I
		// =>
		// C = L^-1 L^-T
		// =>
		// C^-T = L L^T
		// => (C is symmetric)
		// C^-1 = L L^T
		// =>
		// L = Cholesky(C^-1)

		const CholeskyDecompositionD invCholesky(
			inverse(covariance));
		
		// Our samples are row vectors so we
		// multiply by the transpose from the right.

		*signal *= transpose(invCholesky.lower());
	}		

}
