#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"

#include <pastel/math/matrix_tools.h>
#include <pastel/math/cholesky_decomposition_tools.h>

#include <pastel/sys/string_tools.h>

using namespace Tim;

namespace
{

	template <typename NormBijection>
	void testMutualInformationCase(
		const std::string& name,
		const SignalPtr& jointSignal,
		integer kNearest,
		real maxRelativeError,
		const NormBijection& normBijection,
		real correctMi)
	{
		const real mi = mutualInformation(
			jointSignal,
			kNearest,
			maxRelativeError,
			normBijection);

		log() << name << ": " << mi
			<< " (" << mi - correctMi << ")"
			<< logNewLine;
	}

	void testMutualInformation()
	{
		log() << "Mutual information estimates: " << logNewLine;

		const integer samples = 10000;
		const integer kNearest = 10;
		const real maxRelativeError = 0;
		//const EuclideanNormBijection<real> normBijection;
		const InfinityNormBijection<real> normBijection;

		log() << "2d correlated gaussian" << logNewLine;

		{
			const integer dimension = 2;

			for (integer i = 0;i < 10;++i)
			{
				MatrixD covariance(dimension, dimension);

				const real r = (real)i / 10;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecompositionD cholesky(
					covariance);

				const real det = determinant(cholesky);
				const real cond = (1 + r) / (1 - r);

				ENSURE(cholesky.succeeded());

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					jointSignal,
					kNearest,
					maxRelativeError,
					normBijection,
					correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky)));
			}
		}

		log() << "nD correlated gaussian cond-det covariances" << logNewLine;
		{
			const integer dimension = 5;
			for (integer i = 0;i < 10;++i)
			{
				MatrixD covariance(dimension, dimension);

				const real cond = 10 - i;
				const real det = 1 + i;

				//const real cond = 2;
				//const real det = 1 + i;

				/*
				const real cond = 1 + i;
				const real det = 2;
				*/

				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				/*
				log() << "cond = " << conditionManhattan(covariance)
					<< ", det = " << determinant(covariance) << logNewLine;
				*/

				const CholeskyDecompositionD cholesky(
					covariance);

				ENSURE(cholesky.succeeded());

				//std::cout << covariance << std::endl;

				//std::cout << cholesky.lower() << std::endl;

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				testMutualInformationCase(
					"Cor.Gauss. det " + realToString(det) +
					" cond " + realToString(cond),
					jointSignal,
					kNearest,
					maxRelativeError,
					normBijection,
					correlatedGaussianMutualInformation(
					diagonalProduct(covariance), determinant(cholesky)));

				/*
				MatrixD sampleCovariance;
				computeCovariance(jointSignal, sampleCovariance);
				std::cout << sampleCovariance << std::endl;
				*/
			}
		}

		log() << "2D correlated gaussian pairwise" << logNewLine;
		{
			const integer dimension = 2;

			for (integer i = 0;i < 10;++i)
			{
				MatrixD covariance(dimension, dimension);

				const real r = (real)i / 10;

				covariance |= 
					1, r,
					r, 1;

				const CholeskyDecompositionD cholesky(
					covariance);

				const real det = determinant(cholesky);
				const real cond = r;

				ENSURE(cholesky.succeeded());

				const SignalPtr jointSignal = 
					generateCorrelatedGaussian(samples, dimension, cholesky);

				MatrixD pairwiseMi;
				mutualInformationNaive(
					jointSignal,
					100,
					pairwiseMi);
				std::cout << pairwiseMi << std::endl;
			}
		}
	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
	}

	CallFunction run(testAdd);

}