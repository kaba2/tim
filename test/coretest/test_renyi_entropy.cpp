#include "estimation.h"

#include "tim/core/renyi_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/device/timer.h>

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

#include <pastel/math/cholesky_decomposition_tools.h>
#include <pastel/math/matrix_tools.h>

using namespace Tim;

namespace
{

	void testRenyiEntropyCase(
		const std::string& name,
		const SignalPtr& signal,
		real q,
		real correct)
	{
		const integer kNearest = 1;
		const real maxRelativeError = 0;
		//const integer timeWindowRadius = signal->samples() / 10;

		SignalPtr signalSet[] = {signal};

		const real estimate = renyiEntropyLps(
			forwardRange(signalSet), q, 
			maxRelativeError);

		/*
		const real estimate = renyiEntropyKl(
			forwardRange(signalSet), q, 
			maxRelativeError, kNearest, 
			normBijection);
		*/

		/*
		std::vector<real> entropySet;
		entropySet.reserve(signal->samples());
		temporalRenyiEntropyKl(
			forwardRange(signalSet), signal->samples() / 10, std::back_inserter(entropySet), 
			maxRelativeError, kNearest);
		const real estimate = std::accumulate(entropySet.begin(), entropySet.end(), (real)0) /
			entropySet.size();
		*/

		/*
		log() << name << ": " << estimate << ", correct: " 
			<< correct << logNewLine;
		*/
		log() << name << ": " <<
			realToString(estimate, 4) << " ("
			<< realToString(100 * relativeError<real>(estimate, correct), 2) << "%)" << logNewLine;
	}

	void testRenyiEntropy()
	{
		Timer timer;
		timer.setStart();

		log() << "Computing Renyi entropies using Leonenko-Pronzato-Savani estimator..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer dimension = 5;
		const integer samples = 10000;
		const real q = 2;

		testRenyiEntropyCase(
			"Gaussian(0, 1)",
			generateGaussian(samples, dimension),
			q,
			gaussianRenyiEntropy(q, dimension, 1));

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

				testRenyiEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + realToString(conditionManhattan(covariance), 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					q,
					gaussianRenyiEntropy(q, dimension, determinant(cholesky)));
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

				testRenyiEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + realToString(conditionManhattan(covariance), 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					q,
					gaussianRenyiEntropy(q, dimension, determinant(cholesky)));
			}
		}

		testRenyiEntropyCase(
			"Uniform(-1, 1)",
			generateUniform(samples, dimension),
			q,
			uniformRenyiEntropy(std::pow((real)2, (real)dimension)));

		timer.store();
		log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("renyi_entropy", testRenyiEntropy);
	}

	CallFunction run(testAdd);

}
