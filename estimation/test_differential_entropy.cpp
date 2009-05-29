#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

using namespace Tim;

namespace
{

	void testDifferentialEntropy()
	{
		const integer dimension = 10;
		const integer points = 10000;
		const integer kNearest = 1;
		const real maxRelativeError = 0;

		log() << "Estimates of differential entropies:" << logNewLine;

		//EuclideanNormBijection<real> normBijection;
		InfinityNormBijection<real> normBijection;
		//ManhattanNormBijection<real> normBijection;

		{
			SignalPtr signal = generateGaussian(points, dimension);
			const real estimate = differentialEntropy(signal, kNearest, maxRelativeError, 
				normBijection);
			
			log() << "Gaussian: " << estimate << ", correct: " 
				<< gaussianDifferentialEntropy(dimension, 1)
				<< logNewLine;
		}

		{
			SignalPtr signal = generateUniform(points, dimension);
			const real estimate = differentialEntropy(signal, kNearest, maxRelativeError, 
				normBijection);
			
			log() << "Uniform: " << estimate << ", correct: " 
				<< uniformDifferentialEntropy(1)
				<< logNewLine;
		}
	}

	void testAdd()
	{
		timTestList().add("differential_entropy", testDifferentialEntropy);
	}

	CallFunction run(testAdd);

}
