#include "tim/core/signal_tools.h"

#include <pastel/sys/constantiterator.h>

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

	TIMCORE SignalPtr split(
		const SignalPtr& signal,
		integer dimensionBegin,
		integer dimensionEnd)
	{
		ENSURE_OP(dimensionBegin, <=, dimensionEnd);
		ENSURE_OP(dimensionBegin, >=, 0);
		ENSURE_OP(dimensionEnd, <=, signal->dimension());

		const integer dimension = dimensionEnd - dimensionBegin;
		const integer samples = signal->samples();

		SignalPtr splitSignal(new Signal(
			signal->samples(), dimension));

		splitSignal->data() = signal->data()(
			Point2i(0, dimensionBegin), Point2i(samples, dimensionEnd));

		return splitSignal;
	}

	TIMCORE SignalPtr merge(
		const SignalPtr& xSignal,
		const SignalPtr& ySignal,
		integer yLag)
	{
		ENSURE_OP(yLag, >=, 0);

		const SignalPtr signalSet[2] = {xSignal, ySignal};
		const integer lagSet[2] = {0, yLag};

		return Tim::merge(
			forwardRange(signalSet),
			forwardRange(lagSet));
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

		result = (signal->data() - outerProduct(mean, VectorConstant<real, Dynamic>(1, samples))) * 
			transpose(signal->data() - outerProduct(mean, VectorConstant<real, Dynamic>(1, samples)));
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
