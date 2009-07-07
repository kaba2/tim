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
			for (integer j = 0;j < samples;++j)
			{
				for (integer i = 0;i < dimension;++i)
				{
					stream << signal.data()(i, j) << " ";
				}
				stream << std::endl;
			}
		}
		else
		{
			for (integer j = 0;j < samples;++j)
			{
				stream << signal.data()(0, j) << ", ";
			}
		}

		stream << signal.name() << std::endl;

		return stream;
	}

	TIMCORE SignalPtr generateUniform(
		integer samples,
		integer dimension)
	{
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);

		SignalPtr signal = SignalPtr(new Signal(dimension, samples));

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

		SignalPtr signal = SignalPtr(new Signal(dimension, samples));

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
		ENSURE1(dimension > 0, dimension);
		ENSURE1(samples >= 0, samples);
		ENSURE(covarianceCholesky.lower().width() == dimension);
		ENSURE(covarianceCholesky.succeeded());

		SignalPtr correlatedGaussian = generateGaussian(samples, dimension);

		correlatedGaussian->data() = covarianceCholesky.lower() * correlatedGaussian->data();

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

		SignalPtr signal = SignalPtr(new Signal(dimension, samples));

		MatrixD::Iterator iter = signal->data().begin();
		const MatrixD::Iterator iterEnd = signal->data().end();

		while(iter != iterEnd)
		{
			*iter = randomGeneralizedGaussian<real>(shape, scale);
			++iter;
		}

		return signal;
	}

	TIMCORE void slice(
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

		Tim::slice(jointSignal, partition, marginalSet);
	}

	TIMCORE void slice(
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

			marginalSet.push_back(Tim::slice(jointSignal, partition[x], 
				partition[x] + marginalDimension));
		}
	}

	TIMCORE SignalPtr slice(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE2(dimensionBegin <= dimensionEnd, 
			dimensionBegin, dimensionEnd);
		ENSURE1(dimensionBegin >= 0, dimensionBegin);
		ENSURE2(dimensionEnd <= signal->dimension(), dimensionEnd, signal->dimension());

		const integer dimension = dimensionEnd - dimensionBegin;

		return SignalPtr(new Signal(
			dimension, signal->samples(), 
			&signal->data()(dimensionBegin, 0)));
	}

	TIMCORE SignalPtr merge(
		const std::vector<SignalPtr>& signalList)
	{
		if (signalList.empty())
		{
			return SignalPtr();
		}

		const integer samples = signalList.front()->samples();

		const integer signals = signalList.size();
		integer jointDimension = 0;
		for (integer i = 0;i < signals;++i)
		{
			jointDimension += signalList[i]->dimension();
			ENSURE2(signalList[0]->samples() == signalList[i]->samples(), 
				signalList[0]->samples(), signalList[i]->samples());
		}

		SignalPtr jointSignal(new Signal(jointDimension, samples));
		
		integer dimensionOffset = 0;

		for (integer i = 0;i < signals;++i)
		{
			std::copy(
				signalList[i]->data().begin(),
				signalList[i]->data().end(),
				jointSignal->data().rowBegin(dimensionOffset));

			dimensionOffset += signalList[i]->dimension();
		}

		return jointSignal;
	}

	TIMCORE SignalPtr merge(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal)
	{
		std::vector<SignalPtr> signalSet;
		signalSet.reserve(2);
		signalSet.push_back(aSignal);
		signalSet.push_back(bSignal);
		return Tim::merge(signalSet);
	}

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& resultPointSet)
	{
		const integer dimension = signal->dimension();
		const integer n = signal->samples();
		
		std::vector<PointD> pointSet;
		pointSet.reserve(n);

		PointD point(ofDimension(dimension));
		for (integer j = 0;j < n;++j)
		{
			for (integer i = 0;i < dimension;++i)
			{
				point[i] = signal->data()(i, j);
			}

			pointSet.push_back(point);
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

		const VectorD mean = sum(transpose(signal->data())) / samples;

		result = (signal->data() - outerProduct(mean, VectorConstant<Dynamic, real>(1, samples))) * 
			transpose(signal->data() - outerProduct(mean, VectorConstant<Dynamic, real>(1, samples)));
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
		
		signal->data() = transpose(invCholesky.lower()) * signal->data();
	}		

}
