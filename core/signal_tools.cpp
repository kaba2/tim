#include "tim/core/signal_tools.h"

#include <pastel/math/cholesky_decomposition.h>
#include <pastel/math/matrix_tools.h>

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
					stream << signal.data()(j, i) << " ";
				}
				stream << std::endl;
			}
		}
		else
		{
			for (integer j = 0;j < samples;++j)
			{
				stream << signal.data()(j, 0) << ", ";
			}
		}

		stream << signal.name() << std::endl;

		return stream;
	}

	TIMCORE void slice(
		const SignalPtr& jointSignal,
		std::vector<SignalPtr>& marginalSet)
	{
		const integer dimension = jointSignal->dimension();

		SmallSet<integer> partition;
		partition.reserve(dimension + 1);

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
			ENSURE_OP(marginalDimension, >, 0);

			marginalSet.push_back(Tim::slice(jointSignal, partition[x], 
				partition[x] + marginalDimension));
		}
	}

	TIMCORE SignalPtr slice(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signal->dimension());

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = signal->samples();

		SignalPtr sliceSignal(new Signal(
			signal->samples(), dimension));

		sliceSignal->data() = signal->data()(
			Point2i(0, dimensionBegin), Point2i(samples, dimensionEnd));

		return sliceSignal;
	}

	TIMCORE SignalPtr merge(
		const SignalPtr& aSignal,
		const SignalPtr& bSignal,
		integer bLag)
	{
		ENSURE_OP(bLag, >=, 0);

		const SignalPtr signalSet[2] = {aSignal, bSignal};
		const integer lag[2] = {0, bLag};

		return Tim::merge(
			(const SignalPtr*)signalSet, 
			2,
			(const integer*)lag);
	}

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		std::vector<PointD>& pointSet)
	{
		Tim::constructPointSet(
			signal,
			0, signal->samples(),
			0, signal->dimension(),
			pointSet);
	}

	TIMCORE void constructPointSet(
		const SignalPtr& signal,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<PointD>& pointSet)
	{
		ENSURE_OP(sampleBegin, <=, sampleEnd);
		ENSURE_OP(sampleBegin, >= , 0);
		ENSURE_OP(sampleEnd, <=, signal->samples());
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signal->dimension());

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = sampleEnd - sampleBegin;
		
		pointSet.resize(samples);
		std::fill(pointSet.begin(), pointSet.end(), PointD());

		if (dimension == 0 ||
			samples == 0)
		{
			return;
		}

		for (integer i = sampleBegin;i < sampleEnd;++i)
		{
			PointD point(ofDimension(dimension),
				withAliasing(&signal->data()(i, dimensionBegin)));

			pointSet[i - sampleBegin].swap(point);
		}
	}

	TIMCORE void constructPointSet(
		const std::vector<SignalPtr>& ensemble,
		integer sampleBegin,
		integer sampleEnd,
		integer dimensionBegin,
		integer dimensionEnd,
		std::vector<PointD>& pointSet)
	{
		// Check that the trials all have the
		// same dimension and the number of samples.

		const integer trials = ensemble.size();
		if (trials == 0)
		{
			pointSet.clear();
			return;
		}

		const integer signalSamples = ensemble.front()->samples();
		const integer signalDimension = ensemble.front()->dimension();
		for (integer i = 0;i < trials;++i)
		{
			ENSURE_OP(ensemble[i]->samples(), ==, signalSamples);
			ENSURE_OP(ensemble[i]->dimension(), ==, signalDimension);
		}

		ENSURE_OP(sampleBegin, <=, sampleEnd);
		ENSURE_OP(sampleBegin, >= , 0);
		ENSURE_OP(sampleEnd, <=, signalSamples);
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signalDimension);

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = sampleEnd - sampleBegin;

		pointSet.resize(samples * trials);
		std::fill(pointSet.begin(), pointSet.end(), PointD());

		if (dimension == 0 ||
			samples == 0 ||
			trials == 0)
		{
			return;
		}

		integer index = 0;
		for (integer j = 0;j < trials;++j)
		{
			const SignalPtr signal = ensemble[j];

			for (integer i = sampleBegin;i < sampleEnd;++i)
			{
				PointD point(ofDimension(dimension),
					withAliasing(&signal->data()(i, dimensionBegin)));

				pointSet[index].swap(point);
				++index;
			}
		}
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
		
		// The samples are row vectors, so we
		// multiply with the transpose from the right.

		signal->data() *= invCholesky.lower();
	}		

}
