#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal.h"

#include <pastel/math/uniformsampling.h>

namespace Tim
{

	void testDifferentialEntropy()
	{
		const integer dimension = 3;
		const integer points = 10000;
		const integer kNearest = 1;

		
		SignalPtr signal = newSignal(dimension, points);

		for (integer i = 0;i < points;++i)
		{
			(*signal)[i] = DynamicPoint(randomVectorGaussian<Dynamic, real>(dimension));
		}

		log() << "Estimates of differential entropies:" << logNewLine;

		{
			const real estimate = differentialEntropy(signal, kNearest, 0, EuclideanNormBijection<real>());
			
			log() << "Gaussian: " << estimate << ", correct: " 
				<< ((real)dimension / 2) * std::log(2 * constantPi<real>() * constantNeper<real>())
				<< logNewLine;
		}

		for (integer i = 0;i < points;++i)

		{
			(*signal)[i] = DynamicPoint(randomVector<Dynamic, real>(dimension));
		}

		{
			const real estimate = differentialEntropy(signal, kNearest, 0, EuclideanNormBijection<real>());
			
			log() << "Uniform: " << estimate << ", correct: 0" << logNewLine;
		}
	}

	void testMutualInformation()
	{
		SignalPtr aSignal = newSignal(1, 10000);
		SignalPtr bSignal = newSignal(1, 10000);

						
	}

	void testAdd()
	{
		timTestList().add("DifferentialEntropy", testDifferentialEntropy);
		//timTestList().add("MutualInformation", testMutualInformation);
	}

	CallFunction run(testAdd);

}