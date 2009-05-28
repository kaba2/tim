#include "estimation.h"

#include "tim/core/differential_entropy.h"
#include "tim/core/signal_tools.h"

using namespace Tim;

namespace
{

	void testMutualInformation()
	{
		DynamicMatrix correlation(2, 2);
		
		correlation(0, 0) = 1;
		correlation(0, 1) = 0.5;
		correlation(1, 0) = 0.5;
		correlation(1, 1) = 1;
		
		SignalPtr jointSignal = generateCorrelatedGaussian(2, 10000, correlation);
		
		std::vector<SignalPtr> marginalSet;
		splitDimensions(jointSignal, marginalSet);

										
	}

	void testAdd()
	{
		timTestList().add("mutual_information", testMutualInformation);
	}

	CallFunction run(testAdd);

}