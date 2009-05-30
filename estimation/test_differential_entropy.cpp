#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/string_tools.h>

using namespace Tim;

namespace
{

	void testDifferentialEntropyCase(
		const std::string& name,
		const SignalPtr& signal,
		real correct)
	{
		//EuclideanNormBijection<real> normBijection;
		InfinityNormBijection<real> normBijection;
		//ManhattanNormBijection<real> normBijection;

		const integer kNearest = 1;
		const real maxRelativeError = 0;

		const real estimate = differentialEntropy(signal, kNearest, maxRelativeError, 
			normBijection);
		
		log() << name << ": " << estimate << ", correct: " 
			<< correct << logNewLine;
	}

	void testDifferentialEntropy()
	{
		log() << "Estimates of differential entropies:" << logNewLine;

		const integer dimension = 10;
		const integer points = 10000;

		testDifferentialEntropyCase(
			"Gaussian",
			generateGaussian(points, dimension),
			gaussianDifferentialEntropy(dimension, 1));

		testDifferentialEntropyCase(
			"Uniform",
			generateUniform(points, dimension),
			uniformDifferentialEntropy(1));

		for (integer i = 2;i < 10;i += 2)
		{
			const real shape = i;
			const real scale = std::sqrt((real)i);

			testDifferentialEntropyCase(
				"Generalized gaussian " + integerToString(i) + ", sqrt(" + integerToString(i) + ")",
				generateGeneralizedGaussian(points, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}
	}

	void testAdd()
	{
		timTestList().add("differential_entropy", testDifferentialEntropy);
	}

	CallFunction run(testAdd);

}
