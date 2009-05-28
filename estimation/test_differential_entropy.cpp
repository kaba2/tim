#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

using namespace Tim;

namespace
{

	void testDifferentialEntropy()
	{
		const integer dimension = 3;
		const integer points = 10000;
		const integer kNearest = 1;

		log() << "Estimates of differential entropies:" << logNewLine;

		{
			SignalPtr signal = generateGaussian(dimension, points);
			const real estimate = differentialEntropy(signal, kNearest, 0, EuclideanNormBijection<real>());
			
			log() << "Gaussian: " << estimate << ", correct: " 
				<< ((real)dimension / 2) * std::log(2 * constantPi<real>() * constantNeper<real>())
				<< logNewLine;
		}

		{
			SignalPtr signal = generateUniform(dimension, points);
			const real estimate = differentialEntropy(signal, kNearest, 0, EuclideanNormBijection<real>());
			
			log() << "Uniform: " << estimate << ", correct: 0" << logNewLine;
		}
	}

	void testAdd()
	{
		timTestList().add("differential_entropy", testDifferentialEntropy);
	}

	CallFunction run(testAdd);

}
