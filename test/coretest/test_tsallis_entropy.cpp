#include "estimation.h"

#include "tim/core/tsallis_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/random.h>
#include <pastel/sys/string/string_algorithms.h>

#include <pastel/math/cholesky_decomposition.h>
#include <pastel/math/matrix_determinant.h>
#include <pastel/math/random_matrix.h>

using namespace Tim;

namespace
{

	void testTsallisEntropyCase(
		const std::string& name,
		const Signal& signal,
		dreal q,
		dreal correct)
	{
		//const integer kNearest = 1;
		//const integer timeWindowRadius = signal->samples() / 10;

		Signal signalSet[] = {signal};

		const dreal estimate = tsallisEntropyLps(
			range(signalSet), q);

		/*
		const dreal estimate = tsallisEntropyKl(
			range(signalSet), q, 
			kNearest, 
			norm);
		*/

		/*
		std::vector<dreal> entropySet;
		entropySet.reserve(signal->samples());
		temporalTsallisEntropyKl(
			range(signalSet), signal->samples() / 10, std::back_inserter(entropySet), 
			kNearest);
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

	void testTsallisEntropy()
	{
		//Timer timer;
		//timer.setStart();

		log() << "Computing Tsallis entropies using Leonenko-Pronzato-Savani estimator..." << logNewLine;
		log() << "Relative errors to correct analytic results shown in brackets." << logNewLine;

		const integer dimension = 5;
		const integer samples = 10000;
		const dreal q = 2;

		testTsallisEntropyCase(
			"Gaussian(0, 1)",
			generateGaussian(samples, dimension),
			q,
			gaussianTsallisEntropy(q, dimension, 1));

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

				testTsallisEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + 
					realToString(cond, 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					q,
					gaussianTsallisEntropy(q, dimension, determinant(cholesky)));
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

				testTsallisEntropyCase(
					"Cor.G.(" + realToString(determinant(covariance), 2) + ", " + 
					realToString(cond, 2) + ")",
					generateCorrelatedGaussian(samples, dimension, cholesky),
					q,
					gaussianTsallisEntropy(q, dimension, determinant(cholesky)));
			}
		}

		testTsallisEntropyCase(
			"Uniform(-1, 1)",
			generateUniform(samples, dimension),
			q,
			uniformTsallisEntropy(q, std::pow((dreal)2, (dreal)dimension)));

		//timer.store();
		//log() << "Finished in " << timer.seconds() << " seconds." << logNewLine;
	}

	void testAdd()
	{
		timTestList().add("tsallis_entropy", testTsallisEntropy);
	}

	CallFunction run(testAdd);

}
