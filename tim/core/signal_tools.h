// Description: Algorithms for Signal's

#ifndef TIM_SIGNAL_TOOLS_H
#define TIM_SIGNAL_TOOLS_H

#include "tim/core/signal.h"
#include "tim/core/signal_merge.h"
#include "tim/core/signal_split.h"
#include "tim/core/signal_properties.h"

#include <pastel/math/matrix/matrix.h>
#include <pastel/math/matrix/matrix_inverse.h>
#include <pastel/math/matrix/cholesky_decomposition.h>

#include <pastel/sys/range.h>
#include <pastel/sys/array/array.h>
#include <pastel/sys/string/string_algorithms.h>

#include <iostream>

namespace Tim
{

	//! Prints the signal into an output stream.
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

	//! Computes the covariance of the signal samples.
	TIM void computeCovariance(
		const Signal& signal,
		MatrixView<dreal>& result)
	{
		integer dimension = signal.dimension();
		integer samples = signal.samples();

		ENSURE_OP(result.rows(), ==, signal.dimension());
		ENSURE_OP(result.cols(), ==, signal.dimension());

		VectorD mean = asVector(asMatrix(signal.data()).colwise().sum() / samples);

		asMatrix(result) = 
			(signal.matrix() - asColumnMatrix(mean).replicate(1, samples)) * 
			(signal.matrix() - asColumnMatrix(mean).replicate(1, samples)).transpose();
		asMatrix(result) /= samples;
	}

	//! Transforms the given signal to identity covariance.
	TIM void normalizeCovariance(
		Signal& signal,
		const Matrix<dreal>& covariance)
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

		Matrix<dreal> invCovariance = inverse(covariance);

		CholeskyDecompositionInplace<dreal> invCholesky(view(invCovariance));
		
		// The samples are row vectors, so we
		// multiply with the transpose from the right.

		signal.matrix() *= asMatrix(invCholesky.lower());
	}		

}

#include "tim/core/signal_generate.h"

#endif
