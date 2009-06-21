#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"
#include "tim/core/mutual_information.h"

using namespace Tim;

namespace
{

	void testMutualInformation()
	{
		const integer dimension = 2;
		const integer size = 10000;
		const integer kNearest = 1;
		const real maxRelativeError = 0;
		const EuclideanNormBijection<real> normBijection;

		log() << "Mutual information estimates: " << logNewLine;

		{
			const real p = 0.5;
			DynamicMatrix correlation(2, 2);
			correlation(0, 0) = 1;
			correlation(1, 0) = p;
			correlation(0, 1) = p;
			correlation(1, 1) = 1;

			CholeskyDecomposition<Dynamic, real> cholesky(
				correlation);
			
			SignalPtr jointSignal = generateCorrelatedGaussian(size, dimension, cholesky);
			const real mi = mutualInformation(
				jointSignal,
				kNearest,
				maxRelativeError,
				normBijection);
			log() << "Correlated gaussians: " << mi
				<< ", correct: " << correlatedGaussianMutualInformation(correlation)
				<< logNewLine;
		}

	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
	}

	CallFunction run(testAdd);

}