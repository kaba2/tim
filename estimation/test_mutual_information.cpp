#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"

#include <pastel/math/matrix_tools.h>
#include <pastel/math/cholesky_decomposition_tools.h>

using namespace Tim;

namespace
{

	void testMutualInformation()
	{
		const integer dimension = 2;
		const integer size = 1000;
		const integer kNearest = 1;
		const real maxRelativeError = 0;
		const EuclideanNormBijection<real> normBijection;

		log() << "Mutual information estimates: " << logNewLine;

		{
			DynamicMatrix covariance(dimension, dimension);

			/*
			const real cond = 2;
			const real det = 1;
			setRandomSymmetricPositiveDefinite(
				det, cond, covariance);
			*/

			const real r = 0.5;
			covariance(0, 0) = 1;
			covariance(1, 0) = r;
			covariance(0, 1) = r;
			covariance(1, 1) = 1;

			const CholeskyDecomposition<Dynamic, real> cholesky(
				covariance);
			
			const SignalPtr jointSignal = 
				generateCorrelatedGaussian(size, dimension, cholesky);
			const real mi = mutualInformation(
				jointSignal,
				kNearest,
				maxRelativeError,
				normBijection);
			log() << "Correlated gaussians: " << mi
				<< ", correct: " << correlatedGaussianMutualInformation(determinant(cholesky))
				<< logNewLine;
		}

	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
	}

	CallFunction run(testAdd);

}