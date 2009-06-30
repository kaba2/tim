#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/device/timer.h>

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/matrix_tools.h>

using namespace Tim;

namespace
{

	void testDifferentialEntropyCase(
		const std::string& name,
		const SignalPtr& signal,
		real correct)
	{
		//EuclideanNormBijection<real> normBijection;
		//InfinityNormBijection<real> normBijection;
		ManhattanNormBijection<real> normBijection;

		const integer kNearest = 1;
		/*
		const real radiusRatio = 
			std::exp((
			lnVolumeUnitSphereInfinity<real>(signal->width()) - 
			normBijection.lnVolumeUnitSphere(signal->width())) / signal->width());
		const real maxRelativeError = 3 * radiusRatio;
		log() << radiusRatio << logNewLine;
		*/
		const real maxRelativeError = 3;

		const real estimate = differentialEntropy(signal, kNearest, maxRelativeError, 
			normBijection);

		/*
		log() << name << ": " << estimate << ", correct: " 
			<< correct << logNewLine;
		*/
		log() << name << ": " <<
			realToString(estimate, 4) << " ("
			<< realToString(100 * relativeError<real>(estimate, correct), 2) << "%)" << logNewLine;
	}

	void testDifferentialEntropy()
	{
		Timer timer;
		timer.setStart();

		//for (integer iter = 0;iter < 50;++iter)
		{

		log() << "Computing differential entropies using Kozachenko-Leonenko estimator..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer dimension = 10;
		const integer samples = 1000;

		testDifferentialEntropyCase(
			"Gaussian(0, 1)",
			generateGaussian(samples, dimension),
			gaussianDifferentialEntropy(dimension, 1));

		if (dimension > 1)
		{
			for (integer i = 1;i <= 32;i *= 2)
			{
				MatrixD covariance(dimension, dimension);
				const real det = (real)i;
				const real cond = 1;
				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				const CholeskyDecompositionD cholesky(
					covariance);

				testDifferentialEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + realToString(conditionManhattan(covariance), 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					gaussianDifferentialEntropy(dimension, determinant(cholesky)));
			}

			for (integer i = 1;i <= 32;i *= 2)
			{
				MatrixD covariance(dimension, dimension);
				const real det = 1;
				const real cond = (real)i;
				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				const CholeskyDecompositionD cholesky(
					covariance);

				testDifferentialEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + realToString(conditionManhattan(covariance), 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					gaussianDifferentialEntropy(dimension, determinant(cholesky)));
			}
		}

		testDifferentialEntropyCase(
			"Uniform(-1, 1)",
			generateUniform(samples, dimension),
			uniformDifferentialEntropy(std::pow((real)2, (real)dimension)));

		for (integer i = 2;i < 40;i += 8)
		{
			const real shape = i;
			const real scale = varianceToGeneralizedGaussianScale<real>(shape, 1);

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 3) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 2;i < 40;i += 8)
		{
			const real shape = i;
			const real scale = 1;

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 2) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 1;i < 10;++i)
		{
			const real shape = 2 - i * 0.1;
			const real scale = 1;

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 2) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}
		}

		timer.store();
		log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("differential_entropy", testDifferentialEntropy);
	}

	CallFunction run(testAdd);

}
