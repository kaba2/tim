#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/random.h>
#include <pastel/sys/string/string_algorithms.h>

#include <pastel/math/cholesky_decomposition.h>
#include <pastel/math/matrix_determinant.h>
#include <pastel/math/random_matrix.h>

using namespace Tim;

namespace
{

	void testDifferentialEntropyCase(
		const std::string& name,
		const Signal& signal,
		dreal correct)
	{
		Euclidean_Norm<dreal> norm;
		//Maximum_Norm<dreal> norm;
		//Manhattan_Norm<dreal> norm;

		//const integer kNearest = 1;
		//const integer timeWindowRadius = signal->samples() / 10;

		Signal signalSet[] = {signal};

		const dreal estimate = differentialEntropyNk(
			range(signalSet), 
			norm);

		/*
		const dreal estimate = differentialEntropyKl(
			range(signalSet), kNearest, 
			norm);
		*/

		/*
		std::vector<dreal> entropySet;
		entropySet.reserve(signal->samples());
		temporalDifferentialEntropyKl(
			range(signalSet), signal->samples() / 10, std::back_inserter(entropySet), 
			kNearest, norm);
		const dreal estimate = std::accumulate(entropySet.begin(), entropySet.end(), (dreal)0) /
			entropySet.size();
		*/

		/*
		log() << name << ": " << estimate << ", correct: " 
			<< correct << logNewLine;
		*/
		log() << name << ": " <<
			realToString(estimate, 4) << " ("
			<< realToString(100 * relativeError<dreal>(estimate, correct), 2) << "%)" << logNewLine;
	}

	void testDifferentialEntropy()
	{
		//Timer timer;
		//timer.setStart();

		//for (integer iter = 0;iter < 50;++iter)
		{

		log() << "Computing differential entropies using Kozachenko-Leonenko estimator..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer dimension = 5;
		const integer samples = 10000;

		testDifferentialEntropyCase(
			"Gaussian(0, 1)",
			generateGaussian(samples, dimension),
			gaussianDifferentialEntropy(dimension, 1));

		if (dimension > 1)
		{
			for (integer i = 1;i <= 32;i *= 2)
			{
				Matrix<dreal> covariance(dimension, dimension);
				const dreal det = (dreal)i;
				const dreal cond = 1;
				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				const CholeskyDecompositionInplace<dreal> cholesky(
					covariance);

				testDifferentialEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + 
					realToString(cond) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					gaussianDifferentialEntropy(dimension, determinant(cholesky)));
			}

			for (integer i = 1;i <= 32;i *= 2)
			{
				Matrix<dreal> covariance(dimension, dimension);
				const dreal det = 1;
				const dreal cond = (dreal)i;
				setRandomSymmetricPositiveDefinite(
					det, cond, covariance);

				const CholeskyDecompositionInplace<dreal> cholesky(
					covariance);

				testDifferentialEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + 
					realToString(cond) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					gaussianDifferentialEntropy(dimension, determinant(cholesky)));
			}
		}

		testDifferentialEntropyCase(
			"Uniform(-1, 1)",
			generateUniform(samples, dimension),
			uniformDifferentialEntropy(std::pow((dreal)2, (dreal)dimension)));

		for (integer i = 2;i < 40;i += 8)
		{
			const dreal shape = i;
			const dreal scale = varianceToGeneralizedGaussianScale<dreal>(shape, 1);

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 3) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 2;i < 40;i += 8)
		{
			const dreal shape = i;
			const dreal scale = 1;

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 2) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 1;i < 10;++i)
		{
			const dreal shape = 2 - i * 0.1;
			const dreal scale = 1;

			testDifferentialEntropyCase(
				"Gen.G.(" + realToString(shape, 2) + ", " + realToString(scale, 2) + ")",
				generateGeneralizedGaussian(samples, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}
		}

		//timer.store();
		//log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("differential_entropy", testDifferentialEntropy);
	}

	CallFunction run(testAdd);

}
