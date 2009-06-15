#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

#include <pastel/sys/random.h>
#include <pastel/sys/string_tools.h>

using namespace Tim;

namespace
{

	void testDifferentialEntropyCase(
		const std::string& name,
		const SignalPtr& signal,
		real correct)
	{
		EuclideanNormBijection<real> normBijection;
		//InfinityNormBijection<real> normBijection;
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
			"Standard gaussian",
			generateGaussian(points, dimension),
			gaussianDifferentialEntropy(dimension, 1));

		testDifferentialEntropyCase(
			"Uniform [-1, 1]",
			generateUniform(points, dimension),
			uniformDifferentialEntropy(std::pow((real)2, (real)dimension)));

		for (integer i = 2;i < 40;i += 8)
		{
			const real shape = i;
			const real scale = varianceToGeneralizedGaussianScale<real>(shape, 1);

			testDifferentialEntropyCase(
				"Generalized gaussian (shape " + integerToString(shape) + 
				", variance 1)",
				generateGeneralizedGaussian(points, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 2;i < 40;i += 8)
		{
			const real shape = i;
			const real scale = 1;

			testDifferentialEntropyCase(
				"Generalized gaussian (shape " + integerToString(shape) + 
				", scale 1)",
				generateGeneralizedGaussian(points, dimension, shape, scale),
				generalizedGaussianDifferentialEntropy(dimension, shape, scale));
		}

		for (integer i = 1;i < 10;++i)
		{
			const real shape = 2 - i * 0.1;
			const real scale = 1;

			testDifferentialEntropyCase(
				"Generalized gaussian (shape " + integerToString(shape) + 
				", scale 1)",
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
