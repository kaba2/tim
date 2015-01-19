#include "tim/core/signal_tools.h"

#include <pastel/sys/iterator/constant_iterator.h>
#include <pastel/sys/string/string_algorithms.h>

#include <pastel/math/matrix/cholesky_decomposition.h>
#include <pastel/math/matrix/matrix_inverse.h>

namespace Tim
{

	TIM std::ostream& operator<<(
		std::ostream& stream, const Signal& signal)
	{
		integer dimension = signal.dimension();
		integer samples = signal.samples();

		if (samples == 0 || dimension == 0)
		{
			std::cout << "[]" << std::endl;
		}
		else
		{
			std::cout << "[";
			for (integer x = 0; x < dimension;++x)
			{
				if (x > 0)
				{
					std::cout << ";" << std::endl;
				}
				// The beginning time instant is
				// transmitted by padding the first samples
				// of the outputted signal with NaNs.
				for (integer y = 0; y < signal.t();++y)
				{
					if (y > 0)
					{
						std::cout << ", ";
					}
					std::cout << "nan";
				}
				// The NaN-padding is followed by the actual
				// signal.
				for (integer y = 0; y < samples;++y)
				{
					if (y + signal.t() > 0)
					{
						std::cout << ", ";
					}
					std::cout << realToString(signal.data()(y, x));
				}
			}
			std::cout << "]";
		}


		/*
		stream << "[";
		for (integer i = 0;i < dimension;++i)
		{
			if (i > 0)
			{
				stream << ";" << std::endl;
			}
			for (integer j = 0;j < samples;++j)
			{
				if (j > 0)
				{
					stream << ", ";
				}
				stream << signal.data()(j, i);
			}
		}
		stream << "]";
		*/

		return stream;
	}

	TIM void computeCovariance(
		const Signal& signal,
		Matrix<real>& result)
	{
		integer dimension = signal.dimension();
		integer samples = signal.samples();

		result.setSize(dimension, dimension);
		result.set(0);

		VectorD mean = sum(signal.data()) / samples;


		result = (signal.data() - outerProduct(mean, VectorConstant<real, Dynamic>(1, samples))) * 
			transpose(signal.data() - outerProduct(mean, VectorConstant<real, Dynamic>(1, samples)));
		result /= samples;
	}

	TIM void normalizeCovariance(
		Signal& signal,
		const Matrix<real>& covariance)
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

		CholeskyDecomposition<real> invCholesky(
			inverse(covariance));
		
		// The samples are row vectors, so we
		// multiply with the transpose from the right.


		signal.data() *= invCholesky.lower();
	}		

}
