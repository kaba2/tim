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
		const integer dimension = signal.dimension();
		const integer samples = signal.samples();
		if (dimension > 1)
		{
			for (integer i = 0;i < samples;++i)
			{
				for (integer j = 0;j < dimension;++j)
				{
					stream << signal.data()[i][j] << " ";
				}
				stream << std::endl;
			}
		}
		else
		{
			for (integer i = 0;i < samples;++i)
			{
				stream << signal.data()[i][0] << ", ";
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
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);

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
		const CholeskyDecompositionD& correlationCholesky)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);
		ENSURE(correlationCholesky.lower().width() == dimension);
		ENSURE(correlationCholesky.succeeded());

		SignalPtr correlatedGaussian = generateGaussian(samples, dimension);
	
		for (integer i = 0;i < samples;++i)
		{
			correlatedGaussian->data()[i] = 
				correlationCholesky.lower() * correlatedGaussian->data()[i];
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

		MatrixD::Iterator iter = signal->data().begin();
		const MatrixD::Iterator iterEnd = signal->data().end();

		while(iter != iterEnd)
		{
			*iter = randomGeneralizedGaussian<real>(shape, scale);
			++iter;
		}

		return signal;
	}

	TIMCORE void splitMarginal(
		const SignalPtr& jointSignal,
		std::vector<SignalPtr>& marginalSet)
	{
		const integer dimension = jointSignal->dimension();

		SmallSet<integer> partition;
		partition.reserve(dimension);

		for (integer i = 0;i <= dimension;++i)
		{
			partition.insert(i);
		}

		Tim::splitMarginal(
			jointSignal,
			partition,
			marginalSet);
	}

	TIMCORE void splitMarginal(
		const SignalPtr& jointSignal,
		const SmallSet<integer>& partition,
		std::vector<SignalPtr>& marginalSet)
	{
		const integer dimension = jointSignal->dimension();
		const integer samples = jointSignal->samples();
		const integer signals = partition.size() - 1;

		for (integer x = 0;x < signals;++x)
		{
			const integer marginalDimension = 
				partition[x + 1] - partition[x];
			ENSURE1(marginalDimension > 0, marginalDimension);

			const SignalPtr signal = SignalPtr(new Signal(samples, marginalDimension));
			signal->data() = jointSignal->data()(
				Range(0, samples - 1), 
				Range(partition[x], partition[x + 1] - 1));

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

		const integer size = signalList.front()->samples();

		const integer signals = signalList.size();
		integer bigDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			bigDimension += signalList[i]->dimension();
			ENSURE2(signalList[0]->samples() == signalList[i]->samples(), 
				signalList[0]->samples(), signalList[i]->samples());
		}

		SignalPtr bigSignal(new Signal(size, bigDimension));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			copy(signalList[i]->data().constView(),
				subView(bigSignal->data().view(), 
				Rectangle2(0, dimensionOffset, size, 
				dimensionOffset + signalList[i]->dimension())));

			dimensionOffset += signalList[i]->dimension();
		}

		return bigSignal;
	}

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& resultPointSet)
	{
		const integer dimension = signal->dimension();
		const integer n = signal->samples();
		
		std::vector<PointD> pointSet;
		pointSet.reserve(n);

		for (integer i = 0;i < n;++i)
		{
			pointSet.push_back(PointD(signal->data()[i]));
		}

		pointSet.swap(resultPointSet);
	}

	TIMCORE void computeCovariance(
		const SignalPtr& signal,
		MatrixD& result)
	{
		const integer dimension = signal->dimension();
		const integer samples = signal->samples();

		result.setSize(dimension, dimension);
		result.set(0);

		const VectorD mean = sum(signal->data()) / samples;

		for (integer i = 0;i < samples;++i)
		{
			result += outerProduct(signal->data()[i] - mean, signal->data()[i] - mean);
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
		// Problem: find an invertible matrix A
		// by which the signal X transforms
		// into a signal Y = AX 
		// having identity covariance.
		//
		// Solution:
		// 
		// Y Y^T = I
		// =>
		// (AX) (AX)^T = I
		// =>
		// A X X^T A^T = I
		// =>
		// A C A^T = I
		// =>
		// C = A^-1 A^-T
		// =>
		// C^-1 = (A^-1 A^-T)^-1
		// =>
		// C^-1 = A^T A

		// One solution is given by:
		// A^T = Cholesky(C^-1)
		// =>
		// A = Cholesky(C^-1)^T

		const CholeskyDecompositionD invCholesky(
			inverse(covariance));
		
		// Our samples are row vectors so we
		// multiply by the A^T from the right.

		signal->data() *= invCholesky.lower();
	}		

}
